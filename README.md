# 2D-Heat-Conduction-Solver
A C++ program that solves the steady-state 2D heat conduction equation in a four-material composite slab using the Jacobi iterative method with a finite difference approach.

## Overview
This project simulates temperature distribution in a solid domain composed of four different materials, each with a unique thermal conductivity (k). The solver uses a ghost-cell method to apply a combination of three common boundary conditions:

- **Dirichlet**: Fixed temperature on the top and bottom surfaces.
- **Neumann**: Adiabatic (zero heat flux) condition on the left surface.
- **Robin**: Convective heat transfer on the right surface.

The code is structured using an object-oriented approach, with a dedicated HeatConductionSolver class to encapsulate all data and methods related to the simulation. This improves code readability, modularity, and maintainability.

## Core Features
- **Finite Difference Method (FDM)**: Discretizes the heat equation on a 2D grid.
- **Jacobi Method**: An iterative numerical technique used to solve the resulting system of linear equations.
- **Mixed Boundary Conditions**: Handles Dirichlet, Neumann, and Robin conditions simultaneously.
- **Harmonic Mean Conductivity**: Correctly calculates the effective thermal conductivity at material interfaces.
- **VTK Output**: Generates output files (.vtk) that can be visualized with software like ParaView to analyze the temperature field and material properties.

## Visualization
The generated .vtk file can be opened with a scientific visualization tool like ParaView. This allows you to view the 2D temperature contour plot and analyze the thermal behavior of the composite slab.
