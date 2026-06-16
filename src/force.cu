#include "sph_all.h"

__global__
void compute_pressure_soundspeed_factor_kernel(Particle* d_particles, int N, double gamma, double GG1, double inv_dim){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    Particle p_i = d_particles[i];
    p_i.pressure = (gamma - 1.0) * p_i.rho * p_i.u;

    p_i.cs = GG1 * sqrt(p_i.u);
    p_i.factor = 1.0 / (1.0 + inv_dim * p_i.h / p_i.rho * p_i.drho_dh);

    d_particles[i] = p_i;
}

void compute_pressure_soundspeed_factor_gpu(SPHSystem *sph) {
    double GG1 = sqrt(sph->gamma * (sph->gamma - 1.0));
    double inv_dim = 1.0 / (double)sph->dim;

    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    compute_pressure_soundspeed_factor_kernel<<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N, sph->gamma, GG1, inv_dim);
    cudaDeviceSynchronize();

    return;
}

__global__
void reset_accelerations_kernel(Particle* d_particles, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    // write in directly
    d_particles[i].ax = 0.0;
    d_particles[i].ay = 0.0;
    d_particles[i].dudt = 0.0;
}


__global__
void compute_force_kernel(
        Particle* d_particles, 
        int N, 
        double cell_size, 
        int num_cells_x, 
        int num_cells_y, 
        double box_size_x, 
        double box_size_y, 
        int* d_head, 
        int* d_next,
        double alpha, 
        double beta, 
        double epsilon) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    Particle p_i = d_particles[i];

    // current cell
    int cx_i = (int)(p_i.x / cell_size);
    int cy_i = (int)(p_i.y / cell_size);

    int checked_cells[9];
    int num_checked = 0;

    // the surrouding 9 blocks
    for (int d_cy = -1; d_cy <= 1; d_cy++) {
        for (int d_cx = -1; d_cx <= 1; d_cx++) {

            int cx = cx_i + d_cx;
            int cy = cy_i + d_cy;

            // y periodic
            if (cy < 0) cy += num_cells_y;
            else if (cy >= num_cells_y) cy -= num_cells_y;

            // x reflective
            if (cx < 0 || cx >= num_cells_x) continue;

            int cell_index = cx + cy * num_cells_x;

            // avoid recompute
            int already_checked = 0;
            for (int c = 0; c < num_checked; c++) {
                if (checked_cells[c] == cell_index) {
                    already_checked = 1;
                    break;
                }
            }
            if (already_checked) continue;

            checked_cells[num_checked++] = cell_index; 

            int j = d_head[cell_index];

            // visit the linked list
            while (j != -1) {
                Particle p_j = d_particles[j];

                // Y-Periodic Minimum Image Convention
                double dy = p_i.y - p_j.y;
                if (dy > 0.5 * box_size_y) dy -= box_size_y;
                else if (dy < -0.5 * box_size_y) dy += box_size_y;

                // usual particle
                if (i != j) { // no self-force
                    Particle ghost_j = p_j;
                    ghost_j.y  = p_i.y - dy;
                    compute_pairwise_physics_gpu(&p_i, &ghost_j, alpha, beta, epsilon);
                }

                // Left Mirror
                if (p_i.x < cell_size) {
                    Particle ghost_j = p_j;
                    ghost_j.x  = -p_j.x;
                    ghost_j.y  = p_i.y - dy;
                    ghost_j.vx = -p_j.vx; // vx returns
                    compute_pairwise_physics_gpu(&p_i, &ghost_j, alpha, beta, epsilon);
                }

                // Right Mirror
                if (p_i.x > box_size_x - cell_size) {
                    Particle ghost_j = p_j;
                    ghost_j.x  = 2.0 * box_size_x - p_j.x;
                    ghost_j.y  = p_i.y - dy;
                    ghost_j.vx = -p_j.vx; // return by wall
                    compute_pairwise_physics_gpu(&p_i, &ghost_j, alpha, beta, epsilon);
                }

                j = d_next[j]; 
            } // end while (j != -1)
        }
    } // end for d_cx, d_cy

    // write back to the GPU global memory!
    d_particles[i] = p_i;
}

__global__
void self_mirror(
        Particle* d_particles, int N, double box_size_x, 
        double alpha, double beta, double epsilon) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    Particle p_i = d_particles[i];

    // Two self-mirror positions: left-wall mirror and right-wall mirror
    double dxs[2] = {2.0 * p_i.x, 2.0 * (p_i.x - box_size_x)};
    double dvxs[2] = {2.0 * p_i.vx, 2.0 * p_i.vx};

    for (int k = 0; k < 2; k++) {
        // Create a ghost particle representing the particle's own mirror image
        Particle ghost_self = p_i;

        ghost_self.x  = p_i.x - dxs[k];
        ghost_self.vx = p_i.vx - dvxs[k];
        // Y-direction perfectly overlaps with itself; relative distance and velocity are 0
        ghost_self.y  = p_i.y;
        ghost_self.vy = p_i.vy; 

        // Invoke the same core physics function
        compute_pairwise_physics_gpu(&p_i, &ghost_self, alpha, beta, epsilon);
    }
    d_particles[i] = p_i;
}
// Computing force, symmetric 2-layers loop with boundary conditions
void compute_force_xreflective_yperiodic_celllist_gpu(SPHSystem *sph) {

    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    // 1. Initialize acceleration and thermal energy evolution rate for all real particles
    reset_accelerations_kernel<<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N);
    cudaDeviceSynchronize(); // not needed actually 

    compute_force_kernel<<<blocksPerGrid, threadsPerBlock>>>(
            sph->d_particles, 
            sph->N, 
            sph->cell_size, 
            sph->num_cells_x, 
            sph->num_cells_y, 
            sph->box_size_x, 
            sph->box_size_y, 
            sph->d_head, 
            sph->d_next,
            sph->alpha, 
            sph->beta, 
            sph->epsilon);
    cudaDeviceSynchronize();

    // 3. Handle Self-mirror Interactions
    // When a particle gets close to the left or right walls, it interacts with its own mirror image.
    // Each particle calculates this independently, so it can be parallelized with OpenMP.
    self_mirror<<<blocksPerGrid, threadsPerBlock>>>(
            sph->d_particles, 
            sph->N, 
            sph->box_size_x,
            sph->alpha, 
            sph->beta, 
            sph->epsilon);


}
