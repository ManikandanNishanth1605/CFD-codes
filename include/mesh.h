#ifndef MESH_H
#define MESH_H

class Mesh2D{
    public:
        int Nx, Ny;
        double Lx, Ly;
        double dx, dy;

        Mesh2D(int Nx_, int Ny_, double Lx_, double Ly_) {
            Nx = Nx_;
            Ny = Ny_;
            Lx = Lx_;
            Ly = Ly_;
            dx = Lx / (Nx - 1);
            dy = Ly / (Ny - 1);
        }
};

#endif