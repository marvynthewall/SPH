#include "density.h"

void compute_density(SPHSystem2D *sph){
/*
 * OpenMP Parallelization Technique:
 * Adding #pragma omp parallel for will automatically split the loop
 * across multiple CPU cores (e.g., with gcc -fopenmp).
 * Because we are using the "Gather" approach, the loop body only modifies
 * `particles[i].density`. Each core is responsible for a different `i`,
 * so they don't modify the same memory location, thus avoiding Race Conditions.
 */
  for (int i = 0; i < sph->N; i++) {
    double local_rho = 0.0;
    double drhodh = 0.0;
    Particle *p_i = &sph->particles[i]; // indicate i-th particle

    // Placeholder for O(N^2) brute-force implementation.
    // Can be replaced with Neighbor list or Tree search acceleration in the
    // future.
    for (int j = 0; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j]; // indicate j-th particle

      // Calculate distance (2D)
      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;
      double r = sqrt(dx * dx + dy * dy);

      // Gather Approach:
      // Find particles that with distnace less than p_i's h
      if (r <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel_2d(r, p_i->h, &W, &dWdr, &dWdh);

        // Formula (5): rho_i = sum_j ( m_j * W_ij(h_i) )
        local_rho += p_j->mass * W;
        // Formula: drho_i/dhi = sum_j ( m_j * (dW_ij(h_i) / dhi))
        drhodh += p_j->mass * dWdh;
      }
    }

    // Return the calculated density to the particle
    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}
