#pragma once

//Geometry
const int N = 16;
const double xLength = 0.1;
const double yLength = 0.1;
const double dx = xLength / static_cast<double> (N + 1);
const double dy = yLength / static_cast<double> (N + 1);

//Material map, thermal conductivity in W/m*K
const double k1 = 50.0, k2 = 100.0, k3 = 150, k4 = 200.0 ;

//Boundary Conditions (Top & Bottom Temperature can be set in main.cpp)
const double h = 50;	// Convective heat transfer coefficient
const double Tf = 50.0; // Ambient temperature for Robin BC

//Solver settings
const double tol = 1e-06;
const int maxIteration = 10000;





