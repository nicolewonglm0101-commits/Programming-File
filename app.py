import streamlit as st
import ctypes
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# --- DEBUGGING SECTION ---
# Get the folder where THIS python file is sitting
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create the full path to the DLL
lib_name = "cpp_lib.dll"
full_path = os.path.join(script_dir, lib_name)

# Print it to the terminal so we can see it
print(f"\n[DEBUG] Python is running in: {os.getcwd()}")
print(f"[DEBUG] Python script is in:  {script_dir}")
print(f"[DEBUG] Looking for DLL at:   {full_path}")
print(f"[DEBUG] Does DLL exist here?  {os.path.exists(full_path)}\n")
# -------------------------

# Load the library
try:
    if not os.path.exists(full_path):
        st.error(f"CRITICAL ERROR: Python cannot find the file at:\n\n`{full_path}`\n\nPlease check the terminal for details.")
        st.stop()
        
    cpp_lib = ctypes.CDLL(full_path)
    
    # Define input/output types
    cpp_lib.calculate_system_loss.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, 
        ctypes.c_double, ctypes.c_double, ctypes.c_double
    ]
    cpp_lib.calculate_system_loss.restype = ctypes.c_double

except OSError as e:
    st.error(f"DLL Load Failed: {e}")
    st.stop()


# --- UI Setup ---
st.title("Pump & System Curve Analyzer")

# 1. Inputs
st.sidebar.header("System Parameters")
temp_c = st.sidebar.number_input("Temperature (°C)", value=25.0)
L_m = st.sidebar.number_input("Pipe Length (m)", value=13.6)
D_m = st.sidebar.number_input("Pipe Diameter (m)", value=0.062, format="%.4f")
static_head = st.sidebar.number_input("Static Head (m)", value=5.0)

st.sidebar.header("Material & Roughness")
mat_options = {"Galvanized Steel": 0.15e-3, "Commercial Steel": 0.045e-3, "Stainless Steel": 0.002e-3, "Smooth": 0.0, "Custom": -1}
mat_choice = st.sidebar.selectbox("Pipe Material", list(mat_options.keys()))
if mat_choice == "Custom":
    roughness = st.sidebar.number_input("Custom Roughness (m)", value=0.001, format="%.5f")
else:
    roughness = mat_options[mat_choice]

st.sidebar.header("Pump Characteristics")
st.sidebar.write("H = a - bQ²")
pump_a = st.sidebar.number_input("Pump 'a' (Shutoff)", value=50.0)
pump_b = st.sidebar.number_input("Pump 'b'", value=40000.0)

st.sidebar.header("Minor Losses")
K_TOTAL = st.sidebar.number_input("Total K Factor", value=5.0)

st.sidebar.header("Analysis")
Q_max = st.sidebar.number_input("Max Flow (m³/s)", value=0.05, format="%.4f")
Q_required = st.sidebar.number_input("Required Flow (m³/s)", value=0.02, format="%.4f")

# 2. Calculation
if st.button("Run Calculation"):
    flow_rates = np.linspace(0.0001, Q_max, 50)
    system_heads = []
    pump_heads = []
    operating_point = None
    min_diff = float('inf')

    for Q in flow_rates:
        area = 3.14159 * (D_m / 2)**2
        V = Q / area
        
        # Call C++
        loss = cpp_lib.calculate_system_loss(temp_c, V, D_m, L_m, roughness, K_TOTAL)
        
        if loss == -1.0:
            st.error("Error: Temperature out of range (0-100 C).")
            st.stop()
            
        sys_h = static_head + loss
        pump_h = pump_a - (pump_b * (Q**2))
        
        system_heads.append(sys_h)
        pump_heads.append(pump_h)
        
        if abs(sys_h - pump_h) < min_diff:
            min_diff = abs(sys_h - pump_h)
            operating_point = (Q, sys_h)

    # 3. Plotting
    fig, ax = plt.subplots()
    ax.plot(flow_rates, system_heads, label="System Curve")
    ax.plot(flow_rates, pump_heads, label="Pump Curve")
    
    found = min_diff < 1.0 and operating_point[1] > 0
    if found:
        ax.plot(operating_point[0], operating_point[1], 'go', label="Operating Point")
        st.success(f"Operating Point: {operating_point[0]:.4f} m³/s @ {operating_point[1]:.2f} m")
        if operating_point[0] >= Q_required:
            st.info("✅ Pump is ADEQUATE.")
        else:
            st.warning("❌ Pump is INADEQUATE.")
    else:
        st.warning("No intersection found. Increase Max Flow Rate.")
        
    ax.set_xlabel("Flow Rate (m³/s)")
    ax.set_ylabel("Head (m)")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)