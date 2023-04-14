#pragma once
#include <iostream>
#include <cmath>

using namespace std;

// Define the Poisson equation
void poisson1d(int n, double* u, double* f)
{
    double h = 1.0 / (n + 1);
    double h2 = h * h;
    for (int i = 0; i < n; i++) {
        double x = (i + 1) * h;
        u[i] = sin(M_PI * x);
        f[i] = -M_PI * M_PI * sin(M_PI * x);
    }
}

// Compute the residual r = f - Au
void residual(int n, double* r, double* u, double* f)
{
    double h = 1.0 / (n + 1);
    double h2 = h * h;
    for (int i = 1; i < n - 1; i++) {
        r[i] = f[i] - (u[i - 1] - 2 * u[i] + u[i + 1]) / h2;
    }
}

// Full weighting restriction operator
void restriction(int n, double* ru, double* r)
{
    for (int i = 1; i < n / 2 - 1; i++) {
        ru[i] = 0.5 * r[2 * i] + 0.25 * (r[2 * i - 1] + r[2 * i + 1]);
    }
}

// Linear interpolation operator
void interpolation(int n, double* u, double* uu)
{
    for (int i = 0; i < n / 2; i++) {
        uu[2 * i] = u[i];
        uu[2 * i + 1] = 0.5 * (u[i] + u[i + 1]);
    }
}

// Weighted Jacobi iteration
void jacobi(int n, double* u, double* f, int num_iters, double omega)
{
    double h = 1.0 / (n + 1);
    double h2 = h * h;
    for (int k = 0; k < num_iters; k++) {
        for (int i = 1; i < n - 1; i++) {
            u[i] += omega * 0.5 * ((u[i - 1] + u[i + 1]) - h2 * f[i]) - omega * u[i];
        }
    }
}

// V-cycle algorithm
void vcycle(int n, double* u, double* f, int nu1, int nu2, double omega)
{
    if (n <= 2) {
        // Solve exactly on coarsest grid
        u[0] = f[0];
        u[1] = f[1];
        u[2] = f[2];
        return;
    }

    // Smooth initial guess with weighted Jacobi
    jacobi(n, u, f, nu1, omega);

    // Compute residual
    double* r = new double[n];
    residual(n, r, u, f);

    // Restrict residual to coarse grid
    int nc = n / 2;
    double* ru = new double[nc];
    restriction(n, ru, r);

    // Solve exactly on coarse grid
    double* uc = new double[nc];
    vcycle(nc, uc, ru, nu1, nu2, omega);

    // Interpolate coarse grid correction and add to current solution
    double* e = new double[n];
    interpolation(nc, uc, e);
    for (int i = 1; i < n - 1; i++) {
u[i] += e[i];
    }

    // Smooth with weighted Jacobi
    jacobi(n, u, f, nu2, omega);

    // Clean up memory
    delete[] r;
    delete[] ru;
    delete[] uc;
    delete[] e;
}

int main()
{
    // Problem size
    int n = 64;

    // Number of pre- and post-smoothing iterations
    int nu1 = 2;
    int nu2 = 2;

    // Relaxation parameter for weighted Jacobi
    double omega = 0.8;

    // Allocate memory for solution and right-hand side
    double* u = new double[n];
    double* f = new double[n];

    // Initialize problem
    poisson1d(n, u, f);

    // Call V-cycle
    vcycle(n, u, f, nu1, nu2, omega);

    // Print solution
    cout << "Solution: ";
    for (int i = 0; i < n; i++) {
        cout << u[i] << " ";
    }
    cout << endl;

    // Clean up memory
    delete[] u;
    delete[] f;

    return 0;
}
