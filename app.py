import streamlit as st
import ctypes
import os
import sys
import subprocess
import platform
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF
import tempfile

# IMPORTANT: must be the first Streamlit call
st.set_page_config(page_title="Numerical Analysis of Pumped Piping System", layout="wide")

# PATH HELPERS (SCRIPT & PYINSTALLER)

def get_base_dir() -> str:
    """Return base directory for resources (script folder or PyInstaller temp folder)."""
    try:
        return sys._MEIPASS  # type: ignore[attr-defined]
    except AttributeError:
        return os.path.dirname(os.path.abspath(__file__))

def get_resource_path(filename: str) -> str:
    return os.path.join(get_base_dir(), filename)

# LOAD LIBRARY 
def load_cpp_library():
    system_os = platform.system()  # "Windows" or "Linux"
    base_dir = get_base_dir()

    # CASE A: LINUX (Streamlit Cloud)
    if system_os == "Linux":
        lib_name = "cpp_lib.so"
        full_path = os.path.join(base_dir, lib_name)

        # Auto-Compile if missing
        if not os.path.exists(full_path):
            with st.spinner("Compiling C++ code for Linux server..."):
                try:
                    cpp_src = os.path.join(base_dir, "calculations.cpp")
                    subprocess.run(
                        ["g++", "-shared", "-fPIC", "-o", full_path, cpp_src],
                        check=True
                    )
                except subprocess.CalledProcessError:
                    st.error("! Failed to compile C++ code.")
                    st.stop()
                except FileNotFoundError:
                    st.error("! 'g++' not found.")
                    st.stop()

        return ctypes.CDLL(full_path)

    # CASE B: WINDOWS (Your Local Computer)
    else:
        lib_name = "cpp_lib.dll"
        full_path = os.path.join(base_dir, lib_name)

        if not os.path.exists(full_path):
            st.error(
                f" ! Library not found at: {full_path}\n"
                "Please compile 'cpp_lib.dll' in your terminal first!"
            )
            st.stop()

        return ctypes.CDLL(full_path)


# CONFIGURE C++ INTERFACE (ONLY ONCE)
def configure_cpp(cpp_lib):
    # 1. System Head Loss
    cpp_lib.calculate_system_head_loss.argtypes = [ctypes.c_double] * 6
    cpp_lib.calculate_system_head_loss.restype = ctypes.c_double

    # 2. Pump Head
    cpp_lib.calculate_pump_head.argtypes = [ctypes.c_double] * 3
    cpp_lib.calculate_pump_head.restype = ctypes.c_double

    # 3. Power Wrapper (Temp, Q, Head, Eff) -> kW
    cpp_lib.calculate_power_wrapper.argtypes = [ctypes.c_double] * 4
    cpp_lib.calculate_power_wrapper.restype = ctypes.c_double

    # 4. NPSH Wrapper (Temp, Q, D_suc, L_suc, Rough, K_suc, Z) -> meters
    cpp_lib.calculate_npsha_wrapper.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double
    ]
    cpp_lib.calculate_npsha_wrapper.restype = ctypes.c_double

    # 5. Details
    cpp_lib.calculate_details.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.POINTER(ctypes.c_double), ctypes.c_int
    ]
    cpp_lib.calculate_details.restype = ctypes.c_int

    # CSV path setter (optional)
    if hasattr(cpp_lib, "set_water_csv_path"):
        csv_path = get_resource_path("water_properties.csv")
        cpp_lib.set_water_csv_path.argtypes = [ctypes.c_char_p]
        if os.path.exists(csv_path):
            cpp_lib.set_water_csv_path(csv_path.encode("utf-8"))
        else:
            st.warning(f" ! CSV Not Found: {csv_path}")

    return cpp_lib

try:
    cpp_lib = configure_cpp(load_cpp_library())
except Exception as e:
    st.error(f"Critical Error loading library: {e}")
    st.stop()


#  UTILITIES
def flow_regime_from_re(Re):
    if Re < 2300:
        return f"Laminar (Re={Re:.0f})"
    elif Re < 4000:
        return f"Transitional (Re={Re:.0f})"
    else:
        return f"Turbulent (Re={Re:.0f})"

# PDF REPORT GENERATOR
def create_pdf_report(inputs, results, fig):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "Pump & System Analysis Report", ln=True, align="C")
    pdf.ln(10)

    # Inputs Section
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "1. System Configuration", ln=True)
    pdf.set_font("Arial", size=10)
    for k, v in inputs.items():
        pdf.cell(0, 6, f"{k}: {v}", ln=True)
    pdf.ln(5)

    # Results Section
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "2. Operational Results", ln=True)
    pdf.set_font("Arial", size=10)
    for k, v in results.items():
        pdf.cell(0, 6, f"{k}: {v}", ln=True)
    pdf.ln(5)

    # Chart
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "3. Performance Curves", ln=True)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp:
        fig.savefig(tmp.name, dpi=100)
        pdf.image(tmp.name, x=10, w=190)

    data = pdf.output(dest="S")
    if isinstance(data, bytearray):
        data = bytes(data)
    elif isinstance(data, str):
        data = data.encode("latin-1")
    return data


# GUI
st.title("Numerical Analysis of Pumped Piping System")

# Keep operating point between reruns (tabs/buttons)
if "operating_point" not in st.session_state:
    st.session_state.operating_point = None  # (Q, H)

# User Manual
with st.expander("User Manual", expanded=False):
    
    st.markdown("""
### User Manual
1. **System Setup:** Enter pipe geometry, static head, temperature, and pipe material/roughness in the sidebar.  
2. **Fittings:** Open **Select Fittings** to add fittings for minor losses. *(Custom K values supported.)*  
3. **Pump Curve:** Enter coefficients **a** and **b**, plus max flow and required flow.  
4. **Tabs:**
   - **System Analysis:** Generate system + pump curves and find operating point.
   - **NPSH Calculation:** Check cavitation risk based on suction piping.
   - **Energy Cost:** Estimate power and yearly electricity cost.

### Assumptions
- The dimensions remains the same throughout the piping system.
- Inlet and outlet tanks are open to the atmosphere.
- The program analyze one-pipe system only.
- Fluid used is water
- Only a single pump is used 
- To conduct the NPSH Calculation and Power Calculation, operating point must be found first in tab 1. Else, default values of Q Required will be used.
- For the NPSH Calculation, if pump is above the tank, enter negative static suction head.
- For NPSH Required, the value can be taken from the pump manufacturer's performance plot at the desired flow rate. 
- For NPSH Calculation, the term (a-1)V^2/2g is negligible, where a = kinetic correction factor
- Can refer manufacturer's pump performance plot for more parameters such as efficiency. 
""")

# Sidebar Inputs
st.sidebar.header("1. Fluid Properties")
temp_c = st.sidebar.number_input("Temperature (°C)", value = 25.0, min_value = 0.0, max_value = 100.0)

st.sidebar.header("2. Discharge Piping System")
L_m = st.sidebar.number_input("Pipe Length (m)", value = 13.6, min_value = 0.0, step=1.0)
D_m = st.sidebar.number_input("Pipe Diameter (m)", value = 0.062, min_value = 0.0, format="%.4f", step=0.001)
static_head = st.sidebar.number_input("Static Head (m)", value = 5.0)

mat_options = {
    "Galvanized Steel": 0.15e-3,
    "Commercial Steel": 0.045e-3,
    "Stainless Steel": 0.002e-3,
    "Custom": -1,
}
mat_choice = st.sidebar.selectbox("Pipe Material", list(mat_options.keys()))
roughness = (
    st.sidebar.number_input("Roughness (m)", min_value = 0.000, format="%.5f")
    if mat_choice == "Custom"
    else mat_options[mat_choice]
)

st.sidebar.header("3. Pump Curve")
st.sidebar.caption("Equation: $H = a - bQ^2$")
pump_a = st.sidebar.number_input("Shutoff Head 'a' (m)", min_value = 0.0, value = 50.0)
pump_b = st.sidebar.number_input("Coefficient 'b'", min_value = 0.0, value = 40000.0)

st.sidebar.header("4. Discharge Fittings")
fittings = {"Threaded Elbow 90": 0.9, "Gate Valve: fully open": 0.2, "Swing Check Valve": 2.0}
calculated_K = 0.0
with st.sidebar.expander("Select Fittings"):
    for name, k in fittings.items():
        qty = st.number_input(f"{name} (K={k})", 0, step=1, key=name)
        calculated_K += qty * k
manual_K = st.sidebar.number_input("Manual K", 0.0)
K_TOTAL = calculated_K + manual_K
st.sidebar.metric("Total K", f"{K_TOTAL:.2f}")

# Defaults
Q_max_default = 0.05
Q_required_default = 0.02

# TABS
tab1, tab2, tab3 = st.tabs([" System Analysis", "NPSH Calculation", "Energy Cost"])

# TAB 1: MAIN ANALYSIS 
with tab1:
    col_a, col_b = st.columns(2)
    with col_a:
        Q_max = st.number_input("Max Flow (m³/s)", value = Q_max_default, min_value = 0.0, format="%.4f")
    with col_b:
        Q_required = st.number_input("Required Flow (m³/s)", value = Q_required_default, min_value = 0.0, format="%.4f")

    if st.button("Run System Analysis", type="primary"):
        # A. Detailed Calculation at Q_required
        details = (ctypes.c_double * 8)()  # density, visc, Re, f, Major, Minor, Total, V
        success = cpp_lib.calculate_details(temp_c, Q_required, D_m, L_m, roughness, K_TOTAL, details, 8)

        if success:
            Re_val, f_val, major_loss, minor_loss = details[2], details[3], details[4], details[5]
            v_calc = details[7]
            regime = flow_regime_from_re(Re_val)

            st.text("Detailed Hydraulics at Q_required:")
            st.code(
                f"""
Q_req        = {Q_required:.4f} m³/s
Velocity     = {v_calc:.4f} m/s
Re           = {Re_val:.2f}
Flow Type    = {regime}
f (friction) = {f_val:.6f}
Major Losses = {major_loss:.4f} m
Minor Losses = {minor_loss:.4f} m
                """.strip(),
                language="text",
            )
        else:
            st.error("Error calculating details. Check inputs.")

        # B. Generate Curves
        Q_vals = np.linspace(0.0001, Q_max, 200)
        sys_H, pump_H = [], []
        min_diff = float("inf")
        operating_point = None

        for q in Q_vals:
            loss = cpp_lib.calculate_system_head_loss(temp_c, q, D_m, L_m, roughness, K_TOTAL)
            if loss == -1:
                st.stop()
            h_s = static_head + loss
            h_p = cpp_lib.calculate_pump_head(pump_a, pump_b, q)

            sys_H.append(h_s)
            pump_H.append(h_p)

            if abs(h_s - h_p) < min_diff:
                min_diff = abs(h_s - h_p)
                operating_point = (q, h_s)

        # Save operating point for other tabs
        st.session_state.operating_point = operating_point

        # C. Plotting
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(Q_vals, sys_H, label="System Curve", linewidth=2)
        ax.plot(Q_vals, pump_H, label="Pump Curve", linestyle="--")
        ax.axvline(Q_required, color="red", alpha=0.5, linestyle=":", label="Required")

        if operating_point:
            op_Q, op_H = operating_point
            ax.plot(op_Q, op_H, "go", label="Op. Point")
            st.success(f"**Operating Point:** {op_Q:.4f} m³/s @ {op_H:.2f} m")

            # D. PDF Report Button
            pdf_data = create_pdf_report(
                {
                    "Temp": f"{temp_c} C",
                    "L": f"{L_m} m",
                    "D": f"{D_m} m",
                    "K": f"{K_TOTAL:.2f}",
                    "Required Q": f"{Q_required}",
                },
                {
                    "Op Flow": f"{op_Q:.4f} m3/s",
                    "Op Head": f"{op_H:.2f} m",
                    "Status": "ADEQUATE" if op_Q >= Q_required else "INADEQUATE",
                },
                fig,
            )
            st.download_button(" Download PDF Report", pdf_data, "Report.pdf", "application/pdf")
        else:
            st.warning("Pump Undersized (No Intersection)")

        ax.grid(True, alpha=0.3)
        ax.legend()
        st.pyplot(fig)

# TAB 2: NPSH 
with tab2:
    st.markdown("NPSH Analysis")
    st.info("$NPSH_A = \\frac{P_{atm}}{\\rho g} + z - \\frac{P_{vap}}{\\rho g} - H_{loss}$")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Suction Line Config:**")
        D_suc = st.number_input("Suction Pipe Diameter (m)", value=D_m, format="%.4f")
        L_suc = st.number_input("Suction Pipe Length (m)", value=2.0)
        K_suc = st.number_input("Suction Fittings, K", value=1.0)
        Rough_suc = st.number_input("Suction Roughness (m)", value=roughness, format="%.5f")
    with col2:
        st.markdown("Conditions:")
        z_suction = st.number_input("Static Suction Head, z (m)", value=-2.0, help="Negative if pump is ABOVE water")
        npsh_req = st.number_input("NPSH Required (m)", value=3.0)

    if st.button("Check Cavitation Risk"):
        op = st.session_state.operating_point
        calc_Q = op[0] if op else Q_required_default

        npsh_a = cpp_lib.calculate_npsha_wrapper(temp_c, calc_Q, D_suc, L_suc, Rough_suc, K_suc, z_suction)

        if npsh_a != -1:
            margin = npsh_a - npsh_req
            c1, c2 = st.columns(2)
            c1.metric("NPSH Available", f"{npsh_a:.2f} m")
            c2.metric("Safety Margin", f"{margin:.2f} m", delta_color="normal" if margin > 0 else "inverse")

            if margin > 0:
                st.success(" Pump Safe from Cavitation")
            else:
                st.error(" ! CAVITATION RISK DETECTED")
        else:
            st.error("Calculation Error (Check Temp/CSV)")

# --- TAB 3: POWER ---
with tab3:
    st.markdown("Power Consumption")
    st.caption("Formula: $Power = \\frac{\\rho g Q H}{\\eta}$")

    c1, c2 = st.columns(2)
    with c1:
        eff_pump = st.slider("Pump Efficiency", min_value = 0.1, max_value = 1.0, value = 0.75)
        elec_cost = st.number_input("Electricity Cost (RM/kWh)", value = 0.10, min_value = 0.0)
    with c2:
        hours = st.number_input("Hours/Year",value = 1000, min_value = 0)
        pump_rated_power = st.number_input("Pump Rated Power (kW)", value = 5.00, min_value = 0.0)

    if st.button("Calculate Power"):
        op = st.session_state.operating_point
        calc_Q = op[0] if op else Q_required_default
        calc_H = op[1] if op else (pump_a - pump_b * calc_Q**2)

        power_kw = cpp_lib.calculate_power_wrapper(temp_c, calc_Q, calc_H, eff_pump)
        cost_year = power_kw * elec_cost * hours

        m1, m2 = st.columns(2)
        m1.metric("Required Brake Horse Power", f"{power_kw:.2f} kW")
        m2.metric("Yearly Cost", f"RM{cost_year:,.2f}")
        
        if pump_rated_power >= power_kw:
            st.success("Pump power is adequate for the required duty.")
        else: 
            st.error("Pump power is not sufficient. Increase the rated power.")
