#ifndef FIELD_H
#define FIELD_H
#include <vector>

class Field2D{
    public:
        int Nx, Ny; // Size of the field 
        std::vector<double> phi; // Field variable

        Field2D(int Nx_, int Ny_){
            // Constructor
            Nx = Nx_;
            Ny = Ny_;
            phi.resize(Nx * Ny, 0.0); // Flattened array
        }

        double& operator()(int i, int j){ 
            // Returns phi(i, j) like a 2D array
            return phi[i + j * Nx];
        }
};

#endif