import streamlit as st
import ctypes
import os
import sys
import subprocess
import platform
import numpy as np
import matplotlib.pyplot as plt

# --- 1. SMART LIBRARY LOADER (Windows + Linux Support) ---
def load_cpp_library():
    # Detect Operating System
    system_os = platform.system() # Returns "Windows" or "Linux"
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # CASE A: LINUX (Streamlit Cloud)
    if system_os == "Linux":
        lib_name = "cpp_lib.so"
        full_path = os.path.join(script_dir, lib_name)
        
        # Check if we need to compile (Auto-Compile)
        if not os.path.exists(full_path):
            with st.spinner("Compiling C++ code for Linux server..."):
                try:
                    subprocess.run(
                        ["g++", "-shared", "-fPIC", "-o", full_path, "calculations.cpp"], 
                        check=True
                    )
                except subprocess.CalledProcessError:
                    st.error("❌ Failed to compile C++ code. Check 'packages.txt' on GitHub.")
                    st.stop()
                except FileNotFoundError:
                    st.error("❌ 'g++' not found. Did you upload 'packages.txt' with 'build-essential'?")
                    st.stop()
        
        return ctypes.CDLL(full_path)

    # CASE B: WINDOWS (Your Local Computer)
    else:
        lib_name = "cpp_lib.dll"
        full_path = os.path.join(script_dir, lib_name)

        if not os.path.exists(full_path):
            st.error(f"⚠️ Library not found at: {full_path}\nPlease compile 'cpp_lib.dll' in your terminal first!")
            st.stop()
            
        return ctypes.CDLL(full_path)

# --- 2. LOAD & CONFIGURE LIBRARY ---
try:
    cpp_lib = load_cpp_library()
    
    # Define input types: temp, V, D, L, roughness, K_total
    cpp_lib.calculate_system_loss.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, 
        ctypes.c_double, ctypes.c_double, ctypes.c_double
    ]
    # Define return type: Head Loss
    cpp_lib.calculate_system_loss.restype = ctypes.c_double

except Exception as e:
    st.error(f"Critical Error loading library: {e}")
    st.stop()


# --- 3. USER INTERFACE (UI) ---
st.title("Pump & System Curve Analyzer")

# Sidebar: System Parameters
st.sidebar.header("1. System Parameters")
temp_c = st.sidebar.number_input("Temperature (°C)", value=25.0)
L_m = st.sidebar.number_input("Pipe Length (m)", value=13.6)
D_m = st.sidebar.number_input("Pipe Diameter (m)", value=0.062, format="%.4f")
static_head = st.sidebar.number_input("Static Head (m)", value=5.0)

# Sidebar: Material
st.sidebar.markdown("---")
st.sidebar.header("2. Material & Roughness")
mat_options = {
    "Galvanized Steel": 0.15e-3, 
    "Commercial Steel": 0.045e-3, 
    "Stainless Steel": 0.002e-3, 
    "Smooth": 0.0, 
    "Custom": -1
}
mat_choice = st.sidebar.selectbox("Pipe Material", list(mat_options.keys()))

if mat_choice == "Custom":
    roughness = st.sidebar.number_input("Custom Roughness (m)", value=0.001, format="%.5f")
else:
    roughness = mat_options[mat_choice]
    st.sidebar.write(f"Roughness used: {roughness} m")

# Sidebar: Pump Data
st.sidebar.markdown("---")
st.sidebar.header("3. Pump Characteristics")
st.sidebar.markdown("Equation: $H_{pump} = a - bQ^2$")
pump_a = st.sidebar.number_input("Pump 'a' (Shutoff Head)", value=50.0)
pump_b = st.sidebar.number_input("Pump 'b' coefficient", value=40000.0)

# Sidebar: Minor Losses
st.sidebar.markdown("---")
st.sidebar.header("4. Minor Losses")
K_TOTAL = st.sidebar.number_input("Total K Factor (Sum of all fittings)", value=5.0)

# Sidebar: Settings
st.sidebar.markdown("---")
st.sidebar.header("5. Analysis Settings")
Q_max = st.sidebar.number_input("Max Flow Rate (m³/s)", value=0.05, format="%.4f")
Q_required = st.sidebar.number_input("Required Flow (m³/s)", value=0.02, format="%.4f")


# --- 4. CALCULATION & PLOTTING ---
if st.button("Run Calculation"):
    
    # Generate flow rate points
    flow_rates = np.linspace(0.0001, Q_max, 50)
    system_heads = []
    pump_heads = []
    
    operating_point = None
    min_diff = float('inf')

    # Loop through each flow rate
    for Q in flow_rates:
        area = 3.14159 * (D_m / 2)**2
        V = Q / area
        
        # --- CALL C++ FUNCTION ---
        loss = cpp_lib.calculate_system_loss(temp_c, V, D_m, L_m, roughness, K_TOTAL)
        
        # Error handling for invalid temperature
        if loss == -1.0:
            st.error("Error: Temperature is out of range (0-100°C) in water_properties.csv.")
            st.stop()
            
        sys_h = static_head + loss
        pump_h = pump_a - (pump_b * (Q**2))
        
        system_heads.append(sys_h)
        pump_heads.append(pump_h)
        
        # Find intersection
        diff = abs(sys_h - pump_h)
        if diff < min_diff:
            min_diff = diff
            operating_point = (Q, sys_h)

    # --- 5. DISPLAY RESULTS ---
    fig, ax = plt.subplots()
    ax.plot(flow_rates, system_heads, label="System Curve", color="blue")
    ax.plot(flow_rates, pump_heads, label="Pump Curve", color="orange")
    
    # Check if a valid intersection was found (tolerance < 1.0m)
    found = min_diff < 1.0 and operating_point[1] > 0
    
    if found:
        # Plot the dot
        ax.plot(operating_point[0], operating_point[1], 'go', label="Operating Point")
        
        # Show text results
        st.success(f"Operating Point: {operating_point[0]:.4f} m³/s @ {operating_point[1]:.2f} m")
        
        if operating_point[0] >= Q_required:
            st.success("✅ Pump is ADEQUATE for the job.")
        else:
            st.error("❌ Pump is INADEQUATE (Flow too low).")
    else:
        st.warning("No intersection found. Try increasing 'Max Flow Rate'.")
        
    ax.set_xlabel("Flow Rate (m³/s)")
    ax.set_ylabel("Head (m)")
    ax.set_title("System Head vs Pump Head")
    ax.legend()
    ax.grid(True)
    
    st.pyplot(fig)
