#include "density.h"

void compute_density(SPHSystem2D *sph) {
  /*
   * OpenMP Parallelization Technique:
   * Adding #pragma omp parallel for will automatically split the loop
   * across multiple CPU cores (e.g., with gcc -fopenmp).
   * Because we are using the "Gather" approach, the loop body only modifies
   * `particles[i].density`. Each core is responsible for a different `i`,
   * so they don't modify the same memory location, thus avoiding Race
   * Conditions.
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

void compute_density_xreflective_yperiodic(SPHSystem2D *sph) {
  for (int i = 0; i < sph->N; i++) {
    double local_rho = 0.0;
    double drhodh = 0.0;
    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j];

      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;

      // Y-Periodic
      if (dy > 0.5 * sph->box_size_y)
        dy -= sph->box_size_y;
      if (dy < -0.5 * sph->box_size_y)
        dy += sph->box_size_y;
      
      double r = sqrt(dx * dx + dy * dy);

      if (r <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel_2d(r, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      // Left mirror (reflected across x=0)
      double dx_L = p_i->x + p_j->x;
      double r_L = sqrt(dx_L * dx_L + dy * dy);
      if (r_L <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel_2d(r_L, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      // Right mirror (reflected across x=box_size_x)
      double dx_R = p_i->x + p_j->x - 2.0 * sph->box_size_x;
      double r_R = sqrt(dx_R * dx_R + dy * dy);
      if (r_R <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel_2d(r_R, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }
    }

    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}
