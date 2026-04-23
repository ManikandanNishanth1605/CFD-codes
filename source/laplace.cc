#include "laplace.h"

int main(){
    Mesh2D mesh(50, 50, 1.0, 1.0);    // Mesh is initiated
    
    Field2D phi(mesh.Nx, mesh.Ny);    // Field is initiated

    Laplace2D solver(mesh, phi);      // Laplace solver is initiated

    solver.applyBoundaryConditions(); // Enforce boundary conditions
    solver.solve();                   // Solve the Laplace equation

    return 0;
}