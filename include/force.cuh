#pragma once
#include "sph_system.h"

__global__ void compute_pressure_soundspeed_factor_kernel(Particle* d_particles, int N, double gamma, double GG1, double inv_dim);

__global__ void reset_accelerations_kernel(Particle* d_particles, int N);

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
        double epsilon);

__device__ __forceinline__ void compute_pairwise_physics_gpu(
        Particle *p_i,
        const Particle *p_j,
        double alpha,
        double beta,
        double epsilon)
{
    double dx = p_i->x - p_j->x;
    double dy = p_i->y - p_j->y;
    double r = sqrt(dx * dx + dy * dy);

    if (r < 1e-12) return; // avoid error

    double max_h = fmax(p_i->h, p_j->h);
    if (r > max_h) return; // Viscosity & Force only active within max_h

    double dvx = p_i->vx - p_j->vx;
    double dvy = p_i->vy - p_j->vy;

    // get the W_ij(h_i) and W_ij(h_j)
    // ⚠️ 注意：cubic_spline_kernel 也必須是 __device__ 函數喔！
    double W_i, dWdr_i, dWdh_i;
    cubic_spline_kernel(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
    double W_j, dWdr_j, dWdh_j;
    cubic_spline_kernel(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

    // --- 1. Pressure Force ---
    double term_i = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
    double term_j = p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;
    double scalar_force = p_j->mass * (term_i + term_j);

    double ax = -scalar_force * (dx / r);
    double ay = -scalar_force * (dy / r);

    p_i->ax += ax;
    p_i->ay += ay;

    // --- 2. Thermal Energy Evolution (Pressure Part) ---
    double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
    double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r);

    p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_j->mass * inner_product_v_dW_i;

    // --- 3. Artificial Viscosity ---
    double r_dot_v = dx * dvx + dy * dvy;
    double Vax = 0.0, Vay = 0.0;

    if (r_dot_v < 0.0) {
        double h_ij = (p_i->h + p_j->h) / 2.0;
        // 使用傳入的 epsilon
        double mu_ij = h_ij * r_dot_v / (r * r + epsilon * (h_ij * h_ij));

        double c_ij = (p_i->cs + p_j->cs) / 2.0;
        double rho_ij = (p_i->rho + p_j->rho) / 2.0;

        // 使用傳入的 alpha, beta
        double PI_ij = (-alpha * c_ij * mu_ij + beta * mu_ij * mu_ij) / rho_ij;

        double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
        p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

        Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dx/r);
        Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dy/r);
        p_i->ax += Vax;
        p_i->ay += Vay;
    }
}
