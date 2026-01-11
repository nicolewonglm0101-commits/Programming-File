import streamlit as st
import ctypes
import os
import sys
import subprocess
import platform

# --- AUTO-COMPILE & LOAD LIBRARY ---
def load_cpp_library():
    # 1. Detect Operating System
    system_os = platform.system() # Returns "Windows" or "Linux"
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # CASE A: LINUX (Streamlit Cloud)
    if system_os == "Linux":
        lib_name = "cpp_lib.so"
        full_path = os.path.join(script_dir, lib_name)
        
        # We MUST compile it because .dll doesn't work here
        if not os.path.exists(full_path):
            with st.spinner("Compiling C++ code for Linux server..."):
                try:
                    # Run the Linux compile command
                    subprocess.run(
                        ["g++", "-shared", "-fPIC", "-o", full_path, "calculations.cpp"], 
                        check=True
                    )
                except subprocess.CalledProcessError:
                    st.error("❌ Failed to compile C++ code on server. Check packages.txt")
                    st.stop()
                except FileNotFoundError:
                    st.error("❌ The server cannot find 'g++'. Did you upload packages.txt?")
                    st.stop()
                    
        return ctypes.CDLL(full_path)

    # CASE B: WINDOWS (Your Laptop)
    else:
        lib_name = "cpp_lib.dll"
        full_path = os.path.join(script_dir, lib_name)

        if not os.path.exists(full_path):
            st.error(f"⚠️ Library not found: {lib_name}\nPlease compile it on your computer first!")
            st.stop()
            
        return ctypes.CDLL(full_path)

# --- LOAD THE LIBRARY ---
try:
    cpp_lib = load_cpp_library()
    
    # Define argument types
    cpp_lib.calculate_system_loss.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, 
        ctypes.c_double, ctypes.c_double, ctypes.c_double
    ]
    cpp_lib.calculate_system_loss.restype = ctypes.c_double

except Exception as e:
    st.error(f"Critical Error loading library: {e}")
    st.stop()

# --- UI Setup ---
st.title("Pump & System Curve Analyzer")
# ... (Paste the rest of your UI code/sliders here) ...
# ... (Make sure you include the inputs and button logic below this point) ...
