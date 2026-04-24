#include "laplace.h"
#include <iostream>

int main(){

    // Inputs
    int Nx = 100, Ny = 100;    // Number of nodes 
    double Lx = 1.0, Ly = 1.0; // Length of domain
    double tol = 1e-6;         // Tolerance
    int maxIter = 10000;       // Maximum iterations
    double omega = 2.1;        // Relaxation factor; (=1) -> Gauss Seidel, (>1 && <2) -> SOR
    std::string filename = "laplace" + std::to_string(Nx) + "_" + std::to_string(Ny) + ".vtk";
    
    // Mesh and Solver logic
    Mesh2D mesh(Nx, Ny, Lx, Ly);        // Mesh is initiated
    Field2D phi(mesh.Nx, mesh.Ny);      // Field is initiated

    Laplace2D solver(mesh, phi, tol, maxIter, omega); // Laplace solver is initiated
    solver.setBoundaryLeft(1.0);        // Enforce boundary condition
    solver.solve();                     // Solve the Laplace equation
    solver.writeVTK(filename);          // VTK output file is written

    return 0;
}