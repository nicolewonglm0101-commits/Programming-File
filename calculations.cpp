#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 

using namespace std;

const double g = 9.81;
const double PI = 3.14159265358979323846;

double calculate_velocity_internal(double flow_rate, double diameter) {
    if (diameter <= 0.0) return 0.0;
    double radius = diameter / 2.0;
    double area = PI * pow(radius, 2.0);
    return flow_rate / area;
}

// Safe Log10 to prevent crashes on invalid inputs
static double safe_log10(double x) {
    if (x <= 0.0) x = 1e-30;
    return log10(x);
}

// FLUID PROPERTIES TABLE
class WaterTable {
private:
    vector<double> T_list;
    vector<double> rho_list;
    vector<double> mu_list;
    vector<double> pvap_list; 
    bool loaded = false;
    string path = "water_properties.csv";

public:
    WaterTable() = default;

    // To Allow changing the CSV path (called from Python)
    void setPath(const string& newPath) {
        path = newPath; 
        loaded = false;
        T_list.clear(); 
        rho_list.clear(); 
        mu_list.clear(); 
        pvap_list.clear();
    }

    // Loads the CSV file into memory
    bool ensureLoaded() {
        if (loaded) return true;

        ifstream file(path);
        if (!file.is_open()) return false;

        string line;
        // Skip Header Line
        if (!getline(file, line)) { file.close(); return false; }

        while (getline(file, line)) {
            if (line.empty()) continue;
            stringstream ss(line);
            string temp;
            
            // Expected Format: T, Rho, Mu, Pvap
            double T_val, Rho_val, Mu_val, Pvap_val = 0.0;

            if (!getline(ss, temp, ',')) break; 
            try { T_val = stod(temp); } catch(...) { continue; }

            if (!getline(ss, temp, ',')) break; 
            try { Rho_val = stod(temp); } catch(...) { continue; }

            if (!getline(ss, temp, ',')) break; 
            try { Mu_val = stod(temp); } catch(...) { continue; }

            // Optional 4th column (Pvap). If missing, default 0.
            if (getline(ss, temp, ',')) {
                 try { Pvap_val = stod(temp); } catch(...) { Pvap_val = 0.0; }
            }

            T_list.push_back(T_val);
            rho_list.push_back(Rho_val);
            mu_list.push_back(Mu_val);
            pvap_list.push_back(Pvap_val);
        }
        file.close();

        if (T_list.size() < 2) return false;
        loaded = true;
        return true;
    }
    
    // Interpolation when not in range
    bool getProps(double T, double &density, double &viscosity, double &pvap) {
        if (!ensureLoaded()) return false;
        
        // Bounds Check
        if (T < T_list.front() || T > T_list.back()) return false;
        
        // Find interval
        for (size_t i = 1; i < T_list.size(); ++i) {
            if (T >= T_list[i-1] && T <= T_list[i]) {
                // Calculate Ratio (0.0 to 1.0)
                double denom = (T_list[i] - T_list[i-1]);
                if (denom == 0) return false; // Prevent divide by zero
                double r = (T - T_list[i-1]) / denom;
                
                // Interpolate
                density   = rho_list[i-1] + (rho_list[i] - rho_list[i-1]) * r;
                viscosity = mu_list[i-1] + (mu_list[i] - mu_list[i-1]) * r;
                pvap      = pvap_list[i-1] + (pvap_list[i] - pvap_list[i-1]) * r;
                
                return true;
            }
        }
        return false;
    }
};

static WaterTable g_waterTable;

// CALCULATION FUNCTION

double Reynolds_Number(double density, double V, double diameter, double dynamic_viscosity) {
    if (dynamic_viscosity == 0.0) return 0.0;
    return (density * V * diameter) / dynamic_viscosity;
}

// Solves Colebrook-White Equation iteratively
double DarcyFrictionFactor(double Re, double diameter, double roughness) {
    if (Re <= 0.0) return NAN;
    
    // Laminar Flow
    if (Re <= 2300.0) {
        return 64.0 / Re;
    } 
    
    // Turbulent / Transitional Flow
    double ed = roughness / diameter;
    
    // Initial Guess ( using Haaland Equation)
    double inv_sqrt_f = -1.8 * safe_log10(pow(ed / 3.7, 1.11) + 6.9 / Re);
    double f = 1.0 / pow(inv_sqrt_f, 2.0);

    // Iterative Refinement
    double f_old;
    const double tolerance = 1e-8;
    const int max_iter = 50;

    for (int i = 0; i < max_iter; ++i) {
        f_old = f;
        double rhs = -2.0 * safe_log10((ed / 3.7) + (2.51 / (Re * sqrt(f_old))));
        f = 1.0 / pow(rhs, 2.0);

        if (abs(f - f_old) < tolerance) break;
    }
    
    return f;
}

double calculate_major_loss(double f, double L, double D, double V) {
    if (D == 0) return 0;
    double part1 = L / D;
    double part2 = pow(V, 2) / (2 * g);
    return f * part1 * part2;
}

double calculate_minor_loss(double K_total, double V) {
    double part1 = pow(V, 2) / (2 * g);
    return K_total * part1;
}


// NPSH
double NPSH_available(double Pvap, double density, double KL_to_pump, double z, 
                      double L_to_pump, double f, double D, double V) { 
    
    const double Patm = 101325.0; // Standard Atm Pressure (Pa)
    // Calculate Total Head Loss in Suction Line
    // HL = (f * L/D + K) * V^2/2g
    double friction_loss = f * (L_to_pump / D);
    double total_res_coeff = friction_loss + KL_to_pump;
    double velocity_head = pow(V, 2) / (2 * g);
    
    double HL_total_to_pump = total_res_coeff * velocity_head; 
 
    // NPSHa Formula
    // NPSHa = (Patm/rho*g) + StaticHead - (Pvap/rho*g) - Losses
    double term_atm = Patm / (density * g);
    double term_vap = Pvap / (density * g);
    
    double NPSH = term_atm + z - term_vap - HL_total_to_pump; 
 
    return NPSH; 
} 
 
// Power Logic
double PumpPower(double density, double Q, double Head, double efficiency) { 
    if (efficiency <= 0.0) efficiency = 0.01; // Prevent div by zero
    
    // Power ( in Watts) = (rho * g * Q * H) / eff
    double Power = (density * g * Q * Head) / efficiency; 
    return Power; 
}

// EXTERNAL INTERFACE

extern "C" {

    // A. Setup
    void set_water_csv_path(const char* p) {
        if (p) g_waterTable.setPath(string(p));
    }

    // B. Main Head Loss Calculator
    double calculate_system_head_loss(double temperature, double flow_rate, double diameter,
                                      double length, double roughness, double K_total) {
        
        // 1. Get Fluid Properties
        double density, viscosity, Pvap;
        if (!g_waterTable.getProps(temperature, density, viscosity, Pvap)) return -1.0;

        // 2. Geometry
        double V = calculate_velocity_internal(flow_rate, diameter);
        
        // 3. FLuid
        double Re = Reynolds_Number(density, V, diameter, viscosity);
        double f  = DarcyFrictionFactor(Re, diameter, roughness);
        
        // 4. Losses
        double HL_major = calculate_major_loss(f, length, diameter, V);
        double HL_minor = calculate_minor_loss(K_total, V);

        return HL_major + HL_minor;
    }

    // C. Pump Curve Calculator
    double calculate_pump_head(double a, double b, double flow_rate) {
        return a - (b * pow(flow_rate, 2.0));
    }

    // D. NPSH Wrapper 
    double calculate_npsha_wrapper(double temp, double flow_rate, double D_suction, 
                                   double L_suction, double roughness, double K_suction, double z_static) {
        
        double density, viscosity, Pvap;
        if (!g_waterTable.getProps(temp, density, viscosity, Pvap)) return -1.0;

        double V = calculate_velocity_internal(flow_rate, D_suction);
        double Re = Reynolds_Number(density, V, D_suction, viscosity);
        double f = DarcyFrictionFactor(Re, D_suction, roughness);

        // Pass calculated f, V, Pvap, Density 
        return NPSH_available(Pvap, density, K_suction, z_static, L_suction, f, D_suction, V);
    }

    // E. Power Wrapper
    double calculate_power_wrapper(double temp, double Q, double Head, double eff) {
        double density, viscosity, Pvap;
        // If temp is invalid, return 0.0
        if (!g_waterTable.getProps(temp, density, viscosity, Pvap)) return 0.0;
        
        double watts = PumpPower(density, Q, Head, eff);
        return watts / 1000.0; // Return in kW
    }

    // F. Detailed Results 
    // Returns 1 if success, 0 if fail
    int calculate_details(double temperature, double flow_rate, double diameter,
                          double length, double roughness, double K_total,
                          double* out, int out_len) {
        
        double density, viscosity, Pvap;
        if (!g_waterTable.getProps(temperature, density, viscosity, Pvap)) return 0;

        double V = calculate_velocity_internal(flow_rate, diameter);
        double Re = Reynolds_Number(density, V, diameter, viscosity);
        double f  = DarcyFrictionFactor(Re, diameter, roughness);
        
        // Safety check for solver failure
        if (!(isfinite(f)) || f <= 0.0) return 0;

        double HL_major = calculate_major_loss(f, length, diameter, V);
        double HL_minor = calculate_minor_loss(K_total, V);
        
        // Populate Python Buffer (Array)
        if (out_len >= 1) out[0] = density;
        if (out_len >= 2) out[1] = viscosity;
        if (out_len >= 3) out[2] = Re;
        if (out_len >= 4) out[3] = f;
        if (out_len >= 5) out[4] = HL_major;
        if (out_len >= 6) out[5] = HL_minor;
        if (out_len >= 7) out[6] = HL_major + HL_minor;
        if (out_len >= 8) out[7] = V;

        return 1;
    }
}
