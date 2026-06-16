#include "sph_all.h"
#include <chrono>

__global__ void compute_timestep_kernel(
        Particle* d_particles, int N, double cfl,
        double cell_size, int num_cells_x, int num_cells_y,
        double box_size_x, double box_size_y,
        const int* d_head, const int* d_next,
        double* d_dt_min) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    Particle p_i = d_particles[i];
    double vmax_i = p_i.cs;

    // --- 使用 Cell List 尋找鄰居 ---
    int cx_i = (int)(p_i.x / cell_size);
    int cy_i = (int)(p_i.y / cell_size);

    for (int d_cy = -1; d_cy <= 1; d_cy++) {
        for (int d_cx = -1; d_cx <= 1; d_cx++) {
            int cx = cx_i + d_cx;
            int cy = cy_i + d_cy;

            if (cy < 0) cy += num_cells_y;
            else if (cy >= num_cells_y) cy -= num_cells_y;
            if (cx < 0 || cx >= num_cells_x) continue;

            int cell_index = cx + cy * num_cells_x;
            int j = d_head[cell_index];

            while (j != -1) {
                Particle p_j = d_particles[j];

                double dy = p_i.y - p_j.y;
                if (dy > 0.5 * box_size_y) dy -= box_size_y;
                else if (dy < -0.5 * box_size_y) dy += box_size_y;

                // (1) 真實粒子
                if (i != j) {
                    double dx = p_i.x - p_j.x;
                    double r = sqrt(dx*dx + dy*dy);
                    if (r > 1.0e-12 && r <= p_i.h) {
                        double dvx = p_i.vx - p_j.vx;
                        double dvy = p_i.vy - p_j.vy;
                        double wij = (dvx * dx + dvy * dy) / r;
                        double vsig_ij = p_i.cs + p_j.cs - ((wij < 0.0) ? 3.0 * wij : 0.0);
                        vmax_i = fmax(vmax_i, vsig_ij);
                    }
                }

                // (2) 左牆壁鏡像粒子 (防止快撞牆時時間步長太大)
                if (p_i.x < cell_size) {
                    double dx = p_i.x - (-p_j.x);
                    double r = sqrt(dx*dx + dy*dy);
                    if (r > 1.0e-12 && r <= p_i.h) {
                        double dvx = p_i.vx - (-p_j.vx);
                        double dvy = p_i.vy - p_j.vy;
                        double wij = (dvx * dx + dvy * dy) / r;
                        double vsig_ij = p_i.cs + p_j.cs - ((wij < 0.0) ? 3.0 * wij : 0.0);
                        vmax_i = fmax(vmax_i, vsig_ij);
                    }
                }

                // (3) 右牆壁鏡像粒子
                if (p_i.x > box_size_x - cell_size) {
                    double dx = p_i.x - (2.0 * box_size_x - p_j.x);
                    double r = sqrt(dx*dx + dy*dy);
                    if (r > 1.0e-12 && r <= p_i.h) {
                        double dvx = p_i.vx - (-p_j.vx);
                        double dvy = p_i.vy - p_j.vy;
                        double wij = (dvx * dx + dvy * dy) / r;
                        double vsig_ij = p_i.cs + p_j.cs - ((wij < 0.0) ? 3.0 * wij : 0.0);
                        vmax_i = fmax(vmax_i, vsig_ij);
                    }
                }

                j = d_next[j];
            }
        }
    }

    // 計算該粒子的 dt，並用原子操作寫入全域最小值
    if (p_i.h > 0.0 && vmax_i > 0.0) {
        double dt_i = cfl * p_i.h / vmax_i;
        if(dt_i>0.01)dt_i=0.01;
        atomicMin_double(d_dt_min, dt_i);
    }
}


// ==========================================
// 2. Host Wrapper (負責記憶體派發與啟動 Kernel)
// ==========================================
double compute_timestep_signal_velocity_gpu(SPHSystem *sph) {
    // 初始化 CPU 端的最小值
    double h_dt_min = DBL_MAX;

    // 在 GPU 上要一小塊空間來存這個 dt_min
    double *d_dt_min;
    cudaMalloc((void**)&d_dt_min, sizeof(double));
    cudaMemcpy(d_dt_min, &h_dt_min, sizeof(double), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    // 發射 Kernel
    compute_timestep_kernel<<<blocksPerGrid, threadsPerBlock>>>(
            sph->d_particles, sph->N, sph->cfl,
            sph->cell_size, sph->num_cells_x, sph->num_cells_y,
            sph->box_size_x, sph->box_size_y,
            sph->d_head, sph->d_next,
            d_dt_min
            );
    cudaDeviceSynchronize();

    // 把算出來的全域最小值抓回 CPU
    cudaMemcpy(&h_dt_min, d_dt_min, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_dt_min); // 記得釋放記憶體

    // 錯誤檢查
    if (h_dt_min == DBL_MAX || !isfinite(h_dt_min) || h_dt_min <= 0.0) {
        fprintf(stderr, "Error: invalid timestep in compute_timestep_signal_velocity. dt=%e\n", h_dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = h_dt_min;
    return h_dt_min;
}


__global__
void Kick_gpu(Particle* d_particles, int N, double dt){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    Particle p_i = d_particles[i];

    p_i.vx += 0.5 * p_i.ax * dt;
    p_i.vy += 0.5 * p_i.ay * dt;

    p_i.u += 0.5 * p_i.dudt * dt;
    p_i.u = fmax(1e-10, p_i.u);

    // write back to VRAM
    d_particles[i] = p_i;
}

__global__
void Drift_gpu(Particle* d_particles, int N, double dt, double box_x, double box_y){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    Particle p_i = d_particles[i];
    p_i.x += p_i.vx * dt;
    p_i.y += p_i.vy * dt;

    // 1. Load data into local variables (CPU registers) to minimize memory access overhead
    double x     = p_i.x;
    double y     = p_i.y;
    double vx    = p_i.vx;

    // 2. Y-Periodic boundary condition (Branchless via ternary operators)
    y -= (y >= box_y) ? box_y : 0.0;
    y += (y < 0.0)    ? box_y : 0.0;

    // 3. X-Reflective boundary condition (Branchless via condition flags)
    // Handle left boundary reflection
    int hit_left = (x < 0.0);
    x  = hit_left ? -x  : x;
    vx = hit_left ? -vx : vx;

    // Handle right boundary reflection (Bug fixed: changed particle x to box_x)
    int hit_right = (x > box_x);
    x  = hit_right ? (2.0 * box_x - x) : x;
    vx = hit_right ? -vx : vx;

    // 4. Store updated states back to the particle structure in memory once
    p_i.x  = x;
    p_i.y  = y;
    p_i.vx = vx;

    d_particles[i] = p_i;
}

double step_leapfrog_kdk_xreflective_yperiodic_gpu(
        SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *)) {
    double dt = calculate_time_step(sph);

    auto t1 = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 256;
    int blocksPerGrid = (sph->N + threadsPerBlock - 1) / threadsPerBlock;

    // Kick
    Kick_gpu<<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N, dt);
    cudaDeviceSynchronize();

    // Drift
    Drift_gpu<<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N, dt, sph->box_size_x, sph->box_size_y);
    cudaDeviceSynchronize();

    sph->time += dt;

    auto t2 = std::chrono::high_resolution_clock::now();

    // Update hydrodynamic quantities at new position
    update_adaptive_h_gpu(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic_celllist_gpu);
    auto t3 = std::chrono::high_resolution_clock::now();
    // compute_density_xreflective_yperiodic(sph);
    compute_density_xreflective_yperiodic_celllist_gpu(sph);
    auto t4 = std::chrono::high_resolution_clock::now();
    compute_pressure_soundspeed_factor_gpu(sph);
    auto t5 = std::chrono::high_resolution_clock::now();
    compute_forces(sph);
    auto t6 = std::chrono::high_resolution_clock::now();

    // Kick
    Kick_gpu<<<blocksPerGrid, threadsPerBlock>>>(sph->d_particles, sph->N, dt);
    cudaDeviceSynchronize();

    auto t7 = std::chrono::high_resolution_clock::now();

    printf("Kick/Drift: %f s | Adapt H: %f s | Density: %f s | Pressure: %f s | Force: %f s | Kick: %f s\n",
            std::chrono::duration<double>(t2 - t1).count(),
            std::chrono::duration<double>(t3 - t2).count(),
            std::chrono::duration<double>(t4 - t3).count(),
            std::chrono::duration<double>(t5 - t4).count(),
            std::chrono::duration<double>(t6 - t5).count(),
            std::chrono::duration<double>(t7 - t6).count()
          );
    return dt;
}
