#pragma once
#include "sph_system.h"

__global__ void  update_adaptive_h_kernel(Particle* d_particles, int N, double eta, double* d_max_change);


__global__ void compute_density_xreflective_yperiodic_celllist_kernel(
    Particle* d_particles,
    int N,
    double box_size_x,
    double box_size_y,
    double cell_size,
    int num_cells_x,
    int num_cells_y,
    const int* d_head,
    const int* d_next
    );
