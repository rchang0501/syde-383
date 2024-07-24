#include "./bits/stdc++.h"

using namespace std;

const double A1 = 0.32 * 0.26;           // Example value for A1
const double A2 = 4.9514 * pow(10, -5);  // Example value for A2
const double g = 9.81;                   // Acceleration due to gravity in m/s^2

const double hi = 0.08;  // the initial height from free surface to the 2 cm mark

const double alpha = 1;            // alpha assume for 1
const double D = 0.00794;          // Example value for D
const double K = 0.5;              // Example value for K
const double epsilon = 0.0000025;  // Example value for pipe roughness
const double rho = 998;            // Example value for fluid density in kg/m^3
const double mu = 0.001;           // Example value for dynamic viscosity in Pa.s
const double tolerance = 1e-5;     // Tolerance for convergence

void print_time(double t) {
    int minutes = static_cast<int>(t) / 60;
    double seconds = t - (minutes * 60);

    cout << fixed << setprecision(2);
    cout << minutes << " minute(s) and " << seconds << " second(s)" << endl;
}

double calculateRe(double rho, double v2, double D, double mu) {
    return (rho * v2 * D) / mu;
}

double calculateF(double epsilon, double D, double Re) {
    return 0.25 / pow(log10((epsilon / D) / 3.7 + 5.74 / pow(Re, 0.9)), 2);
}

double calculateV2(double v2, double h, double L, bool no_loss) {
    double z1 = h + 0.02 + (L / 150);  // height from ground to free surface

    double Re = calculateRe(rho, v2, D, mu);                       // reynolds number
    double f = Re <= 2300 ? 64 / Re : calculateF(epsilon, D, Re);  // choose laminar or turbulent f based on Re

    double numerator = 2 * g * z1;

    double term1 = alpha * (1 - pow((A2 / A1), 2));
    double term2 = (L * f) / D;
    double term3 = K;

    double denominator = term1 + term2 + term3;

    double res = numerator / denominator;
    if (res < 0) {
        cout << "Error: Trying to take sqrt of negative number" << endl;
    }

    return sqrt(res);
}

void time_to_drain(double L) {
    double v1 = 0;
    double v2 = 0.01;  // Initial guess for v2 - 0.53
    double v2_res = 0;

    const double delta_t = 0.1;
    double total_t = 0;

    double h = hi;  // the current height of the free surface

    while (h >= 0) {
        int max_its = 1000;
        int curr_it = 0;
        while (curr_it < max_its) {
            // calculate with loss
            v2_res = calculateV2(v2, h, L, false);
            double error = fabs(v2 - v2_res);

            // cout << "v2: " << v2 << endl;
            // cout << "v2_res: " << v2_res << endl;
            // cout << "error: " << error << endl;
            // cout << endl;

            if (error <= tolerance) {
                // cout << "converged!" << endl;
                break;
            }
            v2 = v2_res;
        }
        v1 = v2_res * (A2 / A1);

        h -= delta_t * v1;
        total_t += delta_t;
    }

    cout << "Length: " << L << ", Time: ";
    print_time(total_t);
    cout << endl;
}

int main() {
    vector<double> lengths = {0.2, 0.3, 0.4, 0.6};

    for (auto L : lengths) {
        time_to_drain(L);
    }

    return 0;
}
