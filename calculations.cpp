#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;

double g = 9.81;

// Helper function to get properties from CSV
bool getWaterProperties(double T, double &density, double &viscosity){
    ifstream file("water_properties.csv");

    if(!file.is_open()){
        return false;
    }

    string line;
    getline(file, line); // Skip header

    double T1, rho1, mu1;
    double T2, rho2, mu2;
    bool first_read = false;

    while(getline(file, line)){
        stringstream ss(line);
        string temp;

        getline(ss, temp, ',');
        double Tcsv = stod(temp);

        getline(ss, temp, ',');
        double rhocsv = stod(temp);

        getline(ss, temp, ',');
        double mucsv = stod(temp);

        if(!first_read){
            T1 = Tcsv;
            rho1 = rhocsv;
            mu1 = mucsv;
            first_read = true;
            continue;
        }

        T2 = Tcsv;
        rho2 = rhocsv;
        mu2 = mucsv;

        if(T >= T1 && T <= T2){
            double ratio = (T - T1) / (T2 - T1);
            density    = rho1 + (rho2 - rho1) * ratio;
            viscosity = mu1  + (mu2  - mu1 ) * ratio;
            file.close();
            return true;
        }

        T1 = T2;
        rho1 = rho2;
        mu1 = mu2;
    }

    file.close();
    return false;
}

double Reynolds_Number (double density, double V, double diameter, double dynamic_viscosity){
    return (density * V * diameter) / dynamic_viscosity;
}

double DarcyFrictionFactor(double Re, double diameter, double roughness){
    double f;
    if(Re <= 2300){
        f = 64 / Re;
    }
    else{
        f = pow(1/(-1.8*log10((6.9/Re)+pow((roughness/diameter)/3.7, 1.11))), 2);
    }
    return f;
}

double calculate_major_loss(double f, double L, double D, double V){
    return f * (L/D) * (pow(V, 2) / (2*g));
}

double calculate_minor_loss(double K_total, double V){
    return K_total * (pow(V, 2) / (2*g));
}

extern "C" {
    double calculate_system_loss(double temperature, double V, double diameter, double length, double roughness, double K_total) {
        double density, viscosity;

        if(!getWaterProperties(temperature, density, viscosity)){
            return -1.0;
        }

        double Re = Reynolds_Number(density, V, diameter, viscosity);
        double f = DarcyFrictionFactor(Re, diameter, roughness);
        double HL = calculate_major_loss(f, length, diameter, V);
        double Hl = calculate_minor_loss(K_total, V);

        return HL + Hl;
    }
}
