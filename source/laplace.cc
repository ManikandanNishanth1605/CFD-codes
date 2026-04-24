#include "laplace.h"
#include <iostream>

int main(){

    // Inputs
    int Nx = 100, Ny = 100; 
    double Lx = 1.0, Ly = 1.0;
    std::string filename = "laplace" + std::to_string(Nx) + "_" + std::to_string(Ny) + ".vtk";
    
    // Mesh and Solver logic
    Mesh2D mesh(Nx, Ny, Lx, Ly);        // Mesh is initiated
    Field2D phi(mesh.Nx, mesh.Ny);      // Field is initiated

    Laplace2D solver(mesh, phi);        // Laplace solver is initiated
    solver.setBoundaryLeft(1.0);        // Enforce boundary condition
    solver.setBoundaryBottom(1.0);        // Enforce boundary condition
    solver.solve();                     // Solve the Laplace equation
    solver.writeVTK(filename);          // VTK output file is written

    return 0;
}