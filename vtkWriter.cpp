#include "vtkWriter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "solverParameters.h"


void exportToVTK(const std::vector<std::vector<double>>& field, const std::string& filename) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write VTK Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Generated from Composite Slab Solver\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";

    // Use the dimensions of the internal grid (N x N), not the ghost cells.
    // The dimensions of the input field vector are (N+2) x (N+2).
   
    vtkFile << "DIMENSIONS " << N << " " << N << " 1\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dy << " 1\n";

    vtkFile << "POINT_DATA " << N * N << "\n";
    vtkFile << "SCALARS field double\n";
    vtkFile << "LOOKUP_TABLE default\n";

    // Write Data Points 
    // Iterate through the entire field, including ghost cells.
    // The VTK file requires data in a flat, row-major format.
    vtkFile << std::fixed << std::setprecision(6);
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            vtkFile << std::fixed << std::setprecision(3) << field[i][j] << " ";
        }
        vtkFile << "\n";
    }

    vtkFile.close();
    std::cout << "Successfully wrote data to " << filename << std::endl;
}
