#include "./bits/stdc++.h"

using namespace std;

double calculateRe(double rho, double Vavg, double D, double mu) {
    return (rho * Vavg * D) / mu;
}

double calculateF(double epsilon, double D, double Re) {
    return 0.25 / pow(log10((epsilon / D) / 3.7 + 5.74 / pow(Re, 0.9)), 2);
}

double calculateV2(double V2, double A2, double A1, double g, double h, double alpha, double L, double Vavg, double D, double K, double epsilon, double rho, double mu) {
    double V1 = V2 * (A2/A1); 
    Vavg = (V1 + V2) / 2;

    // Calculate Reynolds number
    double Re = calculateRe(rho, Vavg, D, mu);
    double f = calculateF(epsilon, D, Re); // use for turbulent
    // double f = 64 / Re; // use for laminar

    double term1 = pow((A2 / A1), 2) * pow(V2, 2);
    double term2 = (2 * g * h) / alpha;
    double term3 = (f * L * pow(V2, 2)) / (D * alpha);
    double term4 = (K * pow(V2, 2)) / alpha;

    double result = term1 + term2 - term3 - term4;

    // Debugging statements
    cout << "V2: " << V2 << ", Re: " << Re << ", f: " << f << endl;
    cout << "term1: " << term1 << ", term2: " << term2 << ", term3: " << term3 << ", term4: " << term4 << endl;
    cout << "result: " << result << endl;
    cout << "V2: " << V2 << endl;
    cout << "V2_result: " << sqrt(result) << endl;
    cout << "delta: " << sqrt(result) - V2 << endl;
    cout << endl;

    if (result < 0) {
        cout << "Error: Trying to take sqrt of negative number" << endl;
        return -1;  // or some other error indication
    }

    return sqrt(result);
}

int main() {
    double A1 = (32 * 26) * pow(10, -4);  // Example value for A1
    double A2 = 1.98 * pow(10, -4);       // Example value for A2
    double g = 9.81;                      // Acceleration due to gravity in m/s^2
    double h = 0.1;                       // Example value for h
    double alpha = 1;                     // alpha turbulent
    // double alpha = 2;                     // alpha laminar
    double L = 0.2;                       // Example value for L
    double Vavg = 3;                      // Example value for Vavg
    double D = 0.00794;                   // Example value for D
    double K = 0.5;                       // Example value for K
    double epsilon = 0.0025;              // Example value for pipe roughness
    double rho = 998;                     // Example value for fluid density in kg/m^3
    double mu = 0.001;                    // Example value for dynamic viscosity in Pa.s
    double tolerance = 1e-6;              // Tolerance for convergence
    double V2 = 0.515;                        // Initial guess for V2 - 0.53

    double V2_result = 0;
    int max_iterations = 100;  // Maximum number of iterations
    int iteration = 0;

    // while (iteration < max_iterations) {
    //     V2_result = calculateV2(V2, A2, A1, g, h, alpha, L, Vavg, D, K, epsilon, rho, mu);
    //     if (V2_result < 0) {
    //         cout << "Error: Calculation resulted in a negative square root term." << endl;
    //         break;
    //     }
    //     if (fabs(V2 - V2_result) < tolerance) {
    //         break;
    //     }
    //     V2 = (V2_result + V2) / 2;
    //     iteration++;
    // }
    
    double best_v2 = 0;
    double best_v2_res = 0;
    double min_error = 100000;
    while(V2 < 0.52) {
        V2_result = calculateV2(V2, A2, A1, g, h, alpha, L, Vavg, D, K, epsilon, rho, mu);
        if (V2_result < 0) {
            cout << "Error: Calculation resulted in a negative square root term." << endl;
            break;
        }

        double error = fabs(V2 - V2_result);

        if (error < tolerance && !isnan(V2) && !isnan(V2_result)) {
            cout << "yay" << endl;
            break;
        }

        if(error < min_error) {
            min_error = error; 
            best_v2 = V2;
            best_v2_res = V2_result;
        }

        V2 += 0.0000001;
    }

    cout << "best: " << best_v2 << ", " << best_v2_res << endl;

    // if (iteration < max_iterations) {
    //     cout << "Converged to V2 = " << V2_result << " in " << iteration << " iterations." << endl;
    // } else {
    //     cout << "Did not converge within the maximum number of iterations." << endl;
    // }

    cout << "V2_result: " << V2_result << endl;
    double Re_res = calculateRe(rho, best_v2, D, mu);
    cout << Re_res << endl;

    return 0;
}