#include "sph_all.h"


// kernel for update_adaptive_h
// particles and d_max_change should be already passed into GPU memory
__global__
void  update_adaptive_h_kernel(Particle* d_particles, int N, double eta, double* d_max_change){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    Particle p_i = d_particles[i];
    double C = eta * eta * p_i.mass;

    double h = p_i.h;
    double rho = p_i.rho;
    double drhodh = p_i.drho_dh;
    double F = rho * h * h - C;
    double dFdh = drhodh * h * h + 2.0 * rho * h;
    if (fabs(dFdh) < 1.0e-14) return;

    double h_new = h - F / dFdh;
    if (!isfinite(h_new) || h_new <= 0.0) return;

    double change = fabs(h_new - h) / h;

    // CPU: if (change > max_change) {max_change = change;}
    // atomic operation for the max_change writing
    atomicMax_double(d_max_change, change);

    d_particles[i].h = h_new;
}

void update_adaptive_h_gpu(SPHSystem *sph, int max_iter, double tol, double eta,
        void (*compute_density_fn)(SPHSystem *)) {
    if (sph == NULL || sph->particles == NULL) {
        fprintf(stderr, "Error: invalid SPH system in update_adaptive_h.\n");
        exit(EXIT_FAILURE);
    }

    if (eta <= 0.0 || !isfinite(eta)) {
        fprintf(stderr, "Error: invalid eta in update_adaptive_h. eta=%e\n", eta);
        exit(EXIT_FAILURE);
    }

    // only max_change is important
    double *d_max_change = NULL;
    cudaMalloc((void**)&d_max_change, sizeof(double));

    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    for (int iter = 0; iter < max_iter; iter++) {

        // in GPU max_change = 0
        double zero = 0.0;
        cudaMemcpy(d_max_change, &zero, sizeof(double), cudaMemcpyHostToDevice);

        // TODO should be sent in GPU version density fn from integrator 
        compute_density_fn(sph);
        // run GPU kernel function
        update_adaptive_h_kernel <<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N, eta, d_max_change);

        double max_change = 0.0;
        // this contains the "sync", the "blocking" method
        cudaMemcpy(&max_change, d_max_change, sizeof(double), cudaMemcpyDeviceToHost);

        if (max_change < tol)break;
    }
    cudaFree(d_max_change);
}


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
        ){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i>=N) return;
    Particle p_i = d_particles[i];
    double local_rho = 0.0;
    double drhodh = 0.0;

    // cell coord
    int cx_i = (int)(p_i.x / cell_size);
    int cy_i = (int)(p_i.y / cell_size);

    int checked_cells[9];
    int num_checked = 0;

    // find the 9 cells surrounding
    for(int d_cy = -1; d_cy<=1; d_cy++){
        for(int d_cx = -1; d_cx <= 1; d_cx ++){
            // the target cell
            int cx = cx_i + d_cx;
            int cy = cy_i + d_cy;

            // deal with the y periodic
            cy += (cy<0) ? num_cells_y : 0;
            cy -= (cy >= num_cells_y) ? num_cells_y : 0;


            // deal with the X reflective
            // ignore the side wall
            if (cx<0 || cx >= num_cells_x)continue;

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

            // TODO d_head ? d_next?
            int j = d_head[cell_index];
            while (j != -1) {
                Particle p_j = d_particles[j];

                // 基礎相對位置
                double dy = p_i.y - p_j.y;

                // Y-Periodic 距離修正 (Minimum Image)
                if (dy > 0.5 * box_size_y) dy -= box_size_y;
                else if (dy < -0.5 * box_size_y) dy += box_size_y;

                // usual particle
                double dx = p_i.x - p_j.x;
                double r = sqrt(dx * dx + dy * dy);

                if (r <= p_i.h) {
                    double W = 0.0, dWdr = 0.0, dWdh = 0.0;
                    cubic_spline_kernel(r, p_i.h, &W, &dWdr, &dWdh);
                    local_rho += p_j.mass * W;
                    drhodh += p_j.mass * dWdh;
                }

                // left mirror in x
                if (p_i.x < p_i.h) {
                    double dx_L = p_i.x + p_j.x; // x_i - (-x_j)
                    double r_L = sqrt(dx_L * dx_L + dy * dy);
                    if (r_L <= p_i.h) {
                        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
                        cubic_spline_kernel(r_L, p_i.h, &W, &dWdr, &dWdh);
                        local_rho += p_j.mass * W;
                        drhodh += p_j.mass * dWdh;
                    }
                }

                // right mirror in x
                if (p_i.x > box_size_x - p_i.h) {
                    double dx_R = p_i.x + p_j.x - 2.0 * box_size_x;
                    double r_R = sqrt(dx_R * dx_R + dy * dy);
                    if (r_R <= p_i.h) {
                        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
                        cubic_spline_kernel(r_R, p_i.h, &W, &dWdr, &dWdh);
                        local_rho += p_j.mass * W;
                        drhodh += p_j.mass * dWdh;
                    }
                }

                // to next particle, if no, the value is -1
                j = d_next[j];
            } // end while
        }
    }
    // write back to the GPU global memory
    d_particles[i].rho = local_rho;
    d_particles[i].drho_dh = drhodh;
}


void compute_density_xreflective_yperiodic_celllist_gpu(SPHSystem *sph) {

    // recalculate the cell list
    // this is in CPU!!
    // move the particles data back to CPU, this happens once for one step
    // TODO possible slow!!
    copy_particles_D2H(sph);

    build_cell_list(sph);

    // trace the total cells in gpu
    // once updated, adjust the memory for the d_head in GPU memory
    static int current_gpu_total_cells = 0;

    if (current_gpu_total_cells != sph->total_cells) {
        if (sph->d_head) {
            cudaFree(sph->d_head);
        }
        cudaMalloc((void**)&sph->d_head, sph->total_cells * sizeof(int));
        current_gpu_total_cells = sph->total_cells;
    }

    // move the new data back
    size_t head_bytes = sph->total_cells * sizeof(int);
    size_t next_bytes = sph->N * sizeof(int);
    cudaMemcpy(sph->d_head, sph->head, head_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(sph->d_next, sph->next, next_bytes, cudaMemcpyHostToDevice);

    // for GPU 
    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    // memory copy
    compute_density_xreflective_yperiodic_celllist_kernel<<<blocksPerGrid, threadsPerBlock>>>(
            sph->d_particles,
            sph->N,
            sph->box_size_x,
            sph->box_size_y,
            sph->cell_size,
            sph->num_cells_x,
            sph->num_cells_y,
            sph->d_head,
            sph->d_next
            );
    cudaDeviceSynchronize();
}
