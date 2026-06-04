#pragma once
#include "sph_system.h"


__global__ void Drift_gpu(Particle* d_particles, int N, double dt, double box_x, double box_y);
__global__ void Kick_gpu(Particle* d_particles, int N, double dt);

__global__ void compute_timestep_kernel(
        Particle* d_particles, int N, double cfl,
        double cell_size, int num_cells_x, int num_cells_y,
        double box_size_x, double box_size_y,
        const int* d_head, const int* d_next,
        double* d_dt_min);
