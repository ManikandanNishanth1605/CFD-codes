#include <iostream>
#include "matplotlibcpp.h"
#include <vector>

namespace plt = matplotlibcpp;
using namespace std;

int main() {
    int N = 100;
    double L = 1.0;
    double dx = L / (N - 1);

    double rho = 1.0;
    double gamma = 1.0;
    double u = 10.0;
    double F = rho * u;
    double D = gamma / dx;

    vector<double> phi(N, 0.0);
    vector<double> x(N, 0.0);
    vector<double> aP(N, 0.0), aE(N, 0.0), aW(N, 0.0), b(N, 0.0);

    // Boundary conditions
    double phi_L = 0.0;
    double phi_R = 1.0;

    phi[0] = phi_L;
    phi[N-1] = phi_R;

    // Position array
    for (int i = 0; i < N; i++) {
        x[i] = i * dx;
    }

    // Assemble coefficients
    for (int i = 1; i < N-1; i++) {
        aE[i] = D;
        aW[i] = D + F;
        aP[i] = aE[i] + aW[i];
        b[i] = 0.0;
    }

    // TDMA coefficients
    vector<double> P(N, 0.0), Q(N, 0.0);

    // Forward sweep
    for (int i = 1; i < N-1; i++) {
        double denom = aP[i] - aW[i] * P[i-1];
        P[i] = aE[i] / denom;
        Q[i] = (b[i] + aW[i] * Q[i-1]) / denom;
    }

    // Back substitution
    for (int i = N-2; i > 0; i--) {
        phi[i] = P[i] * phi[i+1] + Q[i];
    }

    // Output
    for (int i = 0; i < N; i++) {
        cout << i * dx << " " << phi[i] << endl;
    }

    // Plotting results
    plt::plot(x, phi);
    plt::xlabel("x");
    plt::ylabel("phi");
    plt::title("1D Steady-State Convection-Diffusion");
    plt::save("conv-diff-u10_0.png");

    return 0;
}