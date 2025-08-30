#pragma once
#include <vector>
#include <string>

/*
 * Exports a 2D scalar field to a legacy VTK file format.
 * * This function writes the provided 2D vector of doubles to a .vtk file,
 * which can then be opened and visualized in software like ParaView.
 * The file format is a structured points dataset, suitable for a 2D grid.
 * The function assumes a ghost cell layout, where the internal grid size is
 * (N x N) and the total grid size is ((N+2) x (N+2)).
 * * @param field The 2D vector<vector<double>> containing the scalar data to write.
 * @param filename The name of the output .vtk file (e.g., "temperature.vtk").
 */
void exportToVTK(const std::vector<std::vector<double>>& field, const std::string& filename);
