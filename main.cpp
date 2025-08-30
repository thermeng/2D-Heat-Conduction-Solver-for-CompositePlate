#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>
#include "solverParameters.h"
#include "vtkWriter.h"

using namespace std;

class HeatConductionSolver {
public:
    HeatConductionSolver() {
        //Storing field for calculation
        temperature.assign(N + 2, vector<double>(N + 2, 0.0));
        oldTemperature.assign(N + 2, vector<double>(N + 2, 0.0));
        residualField.assign(N + 2, vector<double>(N + 2, 0.0));
        materialMap.assign(N + 2, vector<int>(N + 2, 0));
        kField.assign(N + 2, vector<double>(N + 2, 0.0));
        k_E.assign(N + 2, vector<double>(N + 2));
        k_W.assign(N + 2, vector<double>(N + 2));
        k_N.assign(N + 2, vector<double>(N + 2));
        k_S.assign(N + 2, vector<double>(N + 2));
        a_E.assign(N + 2, vector<double>(N + 2));
        a_W.assign(N + 2, vector<double>(N + 2));
        a_N.assign(N + 2, vector<double>(N + 2));
        a_S.assign(N + 2, vector<double>(N + 2));
        a_P.assign(N + 2, vector<double>(N + 2));

        // Initialize conductivity map
        conductivityMap = {
            {1, k1},
            {2, k2},
            {3, k3},
            {4, k4}
        };

        // Pre-compute material maps and coefficients
        precompute();
    }

    // Main solver function
    void solve() {
        initialize(temperature, 1.0);

        int iteration = 0;
        double residualSum = 0.0, avgChange = 1.0;

        // Initial BC application before the loop starts
        applyBoundaryConditions();

        // Solver Loop
        while (avgChange > tol && iteration < maxIteration) {
            oldTemperature = temperature;

            // Apply boundary conditions to the ghost cells before each Jacobi sweep
            applyBoundaryConditions();

            // Update interior nodes
            for (int i = 1; i <= N; ++i) {
                for (int j = 1; j <= N; ++j) {
                    temperature[i][j] = (a_E[i][j] * oldTemperature[i][j + 1]
                        + a_W[i][j] * oldTemperature[i][j - 1]
                        + a_N[i][j] * oldTemperature[i - 1][j]
                        + a_S[i][j] * oldTemperature[i + 1][j]) / a_P[i][j];
                }
            }

            // Residual calculation
            residualSum = 0.0;
            for (int i = 1; i <= N; ++i) {
                for (int j = 1; j <= N; ++j) {
                    residualField[i][j] = std::abs(oldTemperature[i][j] - temperature[i][j]);
                    residualSum += residualField[i][j];
                }
            }
            avgChange = residualSum / static_cast<double>(N * N);

            cout << "Iteration: " << iteration << " Residual: " << avgChange << '\n';
            ++iteration;
        }
        printResult();
        exportToVTK(temperature, "temperatureField.vtk");
    }

private:
    vector<vector<double>> temperature, oldTemperature, residualField;
    vector<vector<int>> materialMap;
    map<int, double> conductivityMap;
    vector<vector<double>> kField;
    vector<vector<double>> k_E, k_W, k_N, k_S;
    vector<vector<double>> a_E, a_W, a_N, a_S, a_P;

    void initialize(vector<vector<double>>& field, double value) {
        for (size_t i = 0; i < field.size(); ++i) {
            for (size_t j = 0; j < field[i].size(); ++j) {
                field[i][j] = value;
            }
        }
    }

    void applyBoundaryConditions() {
        // Top (Dirichlet)
        for (size_t j = 0; j < temperature[0].size(); ++j) {
            temperature[0][j] = 100;
        }

        // Bottom (Dirichlet)
        for (size_t j = 0; j < temperature[0].size(); ++j) {
            temperature[temperature.size() - 1][j] = 200;
        }

        // Left (Neumann: adiabatic)
        for (size_t i = 0; i < temperature.size(); ++i) {
            temperature[i][0] = temperature[i][1];
        }

        // Right (Robin: convection)
        for (int i = 1; i <= N; ++i) {
            int mat_boundary = materialMap[i][N];
            double k_boundary = conductivityMap.at(mat_boundary);
            temperature[i][N + 1] = (k_boundary * temperature[i][N] + h * dx * Tf) / (k_boundary + h * dx);
        }
    }

    void precompute() {
        // Initialize material map
        for (int i = 0; i < N + 2; ++i) {
            for (int j = 0; j < N + 2; ++j) {
                if (i <= N / 2 && j <= N / 2) materialMap[i][j] = 1;
                else if (i <= N / 2 && j > N / 2) materialMap[i][j] = 2;
                else if (i > N / 2 && j <= N / 2) materialMap[i][j] = 3;
                else materialMap[i][j] = 4;
            }
        }

        // Initial conductivity field for VTK export
        for (int i = 0; i < N + 2; ++i) {
            for (int j = 0; j < N + 2; ++j) {
                kField[i][j] = conductivityMap.at(materialMap[i][j]);
            }
        }
        exportToVTK(kField, "K_conductivity.vtk");

        // Calculate harmonic mean conductivity
        for (int i = 0; i < N + 2; ++i) {
            for (int j = 0; j < N + 2; ++j) {
                int mat_center = materialMap[i][j];
                double k_center = conductivityMap.at(mat_center);

                if (i > 0) {
                    int mat_N = materialMap[i - 1][j];
                    double k_Neighbor = conductivityMap.at(mat_N);
                    k_N[i][j] = (mat_center != mat_N) ?
                        (2.0 * k_center * k_Neighbor / (k_center + k_Neighbor)) : k_center;
                }
                if (i < N + 1) {
                    int mat_S = materialMap[i + 1][j];
                    double k_SNeighbor = conductivityMap.at(mat_S);
                    k_S[i][j] = (mat_center != mat_S) ?
                        (2.0 * k_center * k_SNeighbor / (k_center + k_SNeighbor)) : k_center;
                }
                if (j < N + 1) {
                    int mat_E = materialMap[i][j + 1];
                    double k_ENeighbor = conductivityMap.at(mat_E);
                    k_E[i][j] = (mat_center != mat_E) ?
                        (2.0 * k_center * k_ENeighbor / (k_center + k_ENeighbor)) : k_center;
                }
                if (j > 0) {
                    int mat_W = materialMap[i][j - 1];
                    double k_WNeighbor = conductivityMap.at(mat_W);
                    k_W[i][j] = (mat_center != mat_W) ?
                        (2.0 * k_center * k_WNeighbor / (k_center + k_WNeighbor)) : k_center;
                }
            }
        }

        // Calculate finite difference coefficients
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                a_E[i][j] = k_E[i][j] * dy / dx;
                a_W[i][j] = k_W[i][j] * dy / dx;
                a_N[i][j] = k_N[i][j] * dx / dy;
                a_S[i][j] = k_S[i][j] * dx / dy;
                a_P[i][j] = a_E[i][j] + a_W[i][j] + a_N[i][j] + a_S[i][j];
            }
        }
    }

    void printResult() {
        cout << "-----------Solution-----------\n";
        for (int i = 1; i < static_cast<int>(temperature.size()) - 1; ++i) {
            for (int j = 1; j < static_cast<int>(temperature[i].size()) - 1; ++j) {
                cout << fixed << setprecision(2) << temperature[i][j] << "\t";
            }
            cout << '\n';
        }
    }
};

int main() {
    HeatConductionSolver solver;
    solver.solve();
}