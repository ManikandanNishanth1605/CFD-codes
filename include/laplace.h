#ifndef LAPLACE_H
#define LAPLACE_H
#include "vtk_writer.h"

#include<iostream>
#include<cmath>

class Laplace2D{
    private:
        Mesh2D& mesh;      // Mesh
        Field2D& phi;      // Field
        double tolerance;  // Tolerance
        int maxIter;       // Maximum iterations
        double omega;      // Relaxation factor
    
    public:
        Laplace2D(Mesh2D& mesh_, Field2D& phi_, double tol_, int maxIter_, double omega_) : mesh(mesh_), phi(phi_), tolerance(tol_), maxIter(maxIter_), omega(omega_) {}

        void setInitialValue(double value){
            // Set initial value of mesh for non-boundary nodes
            for (int i = 1; i < mesh.Nx - 1; i++) {
                for (int j = 1; j < mesh.Ny - 1; j++){
                    phi(i, j) = value; 
                }
            }
        }

        void setBoundaryLeft(double value){
            // Left boundary
            for (int j = 0; j < mesh.Ny; j++) {
                phi(0, j) = value; 
            } 
        }
        void setBoundaryRight(double value){
            // Right boundary
            for (int j = 0; j < mesh.Ny; j++) {
                phi(mesh.Nx - 1, j) = value; 
            } 
        }
        void setBoundaryTop(double value){
            // Top boundary
            for (int i = 0; i < mesh.Nx; i++) {
                phi(i, mesh.Ny - 1) = value; 
            } 
        }
        void setBoundaryBottom(double value){
            // Bottom boundary
            for (int i = 0; i < mesh.Nx; i++) {
                phi(i, 0) = value; 
            } 
        }

        void solve(){
            for (int iter = 0; iter < maxIter; iter++){
                double maxErr = 0.0;
                for (int i = 1; i < mesh.Nx - 1; i++){
                    for (int j = 1; j < mesh.Ny - 1; j++){
                        double old_phi = phi(i,j);
                        double phi_GS = 0.25 * (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1)); // Gauss-Seidel step
                        phi(i, j) = (1 - omega) * old_phi + omega * phi_GS;
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