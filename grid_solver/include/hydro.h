#ifndef HYDRO_H
#define HYDRO_H

#include "grid_system.h"

// Calculate time step based on CFL condition
void calculate_dt(GridSystem *grid);

// Apply boundary conditions (0 = Sod 1D, 1 = KH 2D)
void apply_boundary_conditions(GridSystem *grid, int prob_type);

// Perform one hydrodynamics step using MUSCL-Hancock scheme
void hydro_step(GridSystem *grid);

#endif
