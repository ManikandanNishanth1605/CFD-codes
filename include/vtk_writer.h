#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include "mesh.h"
#include "field.h"

#include <vector>
#include <fstream>
#include <string>

class VTKWriter {
public:
    static void writeStructuredGrid2D(
        const std::string& filename,
        Field2D field,
        Mesh2D mesh
    )   // For structured 2D grids 
    {
        std::ofstream file(filename);

        file << "# vtk DataFile Version 3.0\n"; // File version and identifier
        file << "2D Field\n";                   // Header
        file << "ASCII\n";                      // File format (BINARY also an option)
        
        file << "DATASET STRUCTURED_POINTS\n";  // Geometry/topology, structured points 
        file << "DIMENSIONS " << field.Nx << " " << field.Ny << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << mesh.dx << " " << mesh.dy << " 1\n";

        // Dataset Attributes
        file << "POINT_DATA " << field.Nx * field.Ny << "\n";
        file << "SCALARS field double\n";
        file << "LOOKUP_TABLE default\n";

        for (int i = 0; i < field.Nx; i++) {
            for (int j = 0; j < field.Ny; j++) {
                file << field(i, j) << "\n";
            }
        }

        file.close();
    }
};

#endif