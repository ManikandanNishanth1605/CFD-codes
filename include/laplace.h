#ifndef LAPLACE_H
#define LAPLACE_H
#include "vtk_writer.h"

#include<iostream>
#include<cmath>

class Laplace2D{
    private:
        Mesh2D& mesh;
        Field2D& phi;
        double tolerance;
        int maxIter;
    
    public:
        Laplace2D(Mesh2D& mesh_, Field2D& phi_) : mesh(mesh_), phi(phi_), tolerance(1e-6), maxIter(10000) {}

        void applyBoundaryConditions(){
            for (int j = 0; j < mesh.Ny; j++) {
                phi(0, j) = 1.0; // Left boundary
                phi(mesh.Nx - 1, j) = 0.0; // Right boundary
            }
            for (int i = 0; i < mesh.Nx; i++) {
                phi(i, 0) = 0.0; // Bottom boundary
                phi(i, mesh.Ny - 1) = 0.0; // Top boundary
            }
        }

        void solve(){
            for (int iter = 0; iter < maxIter; iter++){
                double maxErr = 0.0;
                for (int i = 1; i < mesh.Nx - 1; i++){
                    for (int j = 1; j < mesh.Ny - 1; j++){
                        double old_phi = phi(i,j);
                        phi(i,j) = 0.25*(phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1));
                        double err = fabs(phi(i,j) - old_phi);

                        if (err > maxErr) {maxErr = err;}
                    }
                }

                if (iter % 100 == 0) {
                    std::cout << "Iteration: " << iter << ", Max Error: " << maxErr << std::endl;
                }

                if (maxErr < tolerance) {
                    std::cout << "Solution converged in " << iter << " iterations" << std::endl;
                    break;               
                }
            }
        }

        void writeVTK(const std::string& filename) {
            VTKWriter::writeStructuredGrid2D(filename, phi, mesh);
        }
        
};

#endif