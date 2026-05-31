#include "density.h"

void update_adaptive_h(SPHSystem *sph, int max_iter, double tol, double eta,
                       void (*compute_density_fn)(SPHSystem *)) {
  if (sph == NULL || sph->particles == NULL) {
    fprintf(stderr, "Error: invalid SPH system in update_adaptive_h.\n");
    exit(EXIT_FAILURE);
  }

  if (eta <= 0.0 || !isfinite(eta)) {
    fprintf(stderr, "Error: invalid eta in update_adaptive_h. eta=%e\n", eta);
    exit(EXIT_FAILURE);
  }

  for (int iter = 0; iter < max_iter; iter++) {
    compute_density_fn(sph);

    double max_change = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(max : max_change) schedule(static)
#endif
    /*Newton iteration: Repeatedly update h until the maximum relative change of
     * h is small enough.*/
    for (int i = 0; i < sph->N; i++) {

      Particle *p_i = &sph->particles[i];
      /*
       * 2D adaptive h target:
       *
       *     rho_i * h_i^2 = eta^2 * m_i
       *
       * Define:
       *
       *     C = eta^2 * m_i
       */
      double C = eta * eta * p_i->mass;

      double h = p_i->h;
      double rho = p_i->rho;
      double drhodh = p_i->drho_dh;
      /*
       * Newton method solves:
       *
       *     F(h) = rho(h) * h^2 - C = 0
       *
       * Its derivative is:
       *
       *     dF/dh = (drho/dh) * h^2 + 2 * rho * h
       */
      double F = rho * h * h - C;
      double dFdh = drhodh * h * h + 2.0 * rho * h;

      /* Avoid division by a very small number.*/
      if (fabs(dFdh) < 1.0e-14) {
        continue;
      }

      double h_new = h - F / dFdh;

      /*Reject invalid smoothing length.*/
      if (!isfinite(h_new) || h_new <= 0.0) {
        continue;
      }

      /*Compute relative change of h for convergence check.*/
      double change = fabs(h_new - h) / h;

      if (change > max_change) {
        max_change = change;
      }

      p_i->h = h_new;
    }
    /*
     * If the largest relative change is smaller than tol,
     * the adaptive h iteration is considered converged.
     */
    if (max_change < tol) {
      printf("adaptive h 2D converged at iter: %d, eta = %.3f\n", iter, eta);
      break;
    }
  }
}

void update_adaptive_h_3d(SPHSystem *sph, int max_iter, double tol, double eta,
                          void (*compute_density_fn)(SPHSystem *)) {
  if (sph == NULL || sph->particles == NULL) {
    fprintf(stderr, "Error: invalid SPH system in update_adaptive_h_3d.\n");
    exit(EXIT_FAILURE);
  }

  if (eta <= 0.0 || !isfinite(eta)) {
    fprintf(stderr, "Error: invalid eta in update_adaptive_h_3d. eta=%e\n",
            eta);
    exit(EXIT_FAILURE);
  }

  for (int iter = 0; iter < max_iter; iter++) {

    compute_density_fn(sph);

    double max_change = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(max : max_change) schedule(static)
#endif
    /*Newton iteration: Repeatedly update h until the maximum relative change of
     * h is small enough.*/
    for (int i = 0; i < sph->N; i++) {

      Particle *p_i = &sph->particles[i];
      /*
       * 3D adaptive h target:
       *
       *     rho_i * h_i^3 = eta^3 * m_i
       *
       * Define:
       *
       *     C = eta^3 * m_i
       */
      double h = p_i->h;
      double rho = p_i->rho;
      double drhodh = p_i->drho_dh;

      /*
       * Newton method solves:
       *
       *     F(h) = rho(h) * h^3 - C = 0
       *
       * Its derivative is:
       *
       *     dF/dh = (drho/dh) * h^3 + 3 * rho * h * h
       */
      double C = eta * eta * eta * p_i->mass;
      double F = rho * h * h * h - C;
      double dFdh = drhodh * h * h * h + 3.0 * rho * h * h;

      /* Avoid division by a very small number.*/
      if (fabs(dFdh) < 1.0e-14) {
        continue;
      }

      double h_new = h - F / dFdh;

      /*Reject invalid smoothing length.*/
      if (!isfinite(h_new) || h_new <= 0.0) {
        continue;
      }

      /*Compute relative change of h for convergence check.*/
      double change = fabs(h_new - h) / h;
      if (change > max_change) {
        max_change = change;
      }

      p_i->h = h_new;
    }
    /*
     * If the largest relative change is smaller than tol,
     * the adaptive h iteration is considered converged.
     */
    if (max_change < tol) {
      printf("adaptive h 3D converged at iter: %d, eta = %.3f\n", iter, eta);
      break;
    }
  }
}

void check_adaptive_h(SPHSystem *sph, double eta, double tol) {
  if (sph == NULL || sph->particles == NULL) {
    fprintf(stderr, "Error: invalid SPH system in check_adaptive_h.\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < sph->N; i++) {

    Particle *p_i = &sph->particles[i];

    double C = eta * eta * p_i->mass;
    double residual = fabs(p_i->rho * p_i->h * p_i->h - C) / C;

    printf("i=%d h=%e rho=%e rho*h^2=%e C=%e residual=%e\n", i, p_i->h,
           p_i->rho, p_i->rho * p_i->h * p_i->h, C, residual);

    if (residual > tol) {
      fprintf(stderr, "Warning: particle %d not converged, residual=%e\n", i,
              residual);
    }
  }
}

void compute_density(SPHSystem *sph) {
/*
 * OpenMP Parallelization Technique:
 * Adding #pragma omp parallel for will automatically split the loop
 * across multiple CPU cores (e.g., with gcc -fopenmp).
 * Because we are using the "Gather" approach, the loop body only modifies
 * `particles[i].rho`. Each core is responsible for a different `i`,
 * so they don't modify the same memory location, thus avoiding Race Conditions.
 */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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
        cubic_spline_kernel(r, p_i->h, &W, &dWdr, &dWdh);

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

void compute_density_3d(SPHSystem *sph) {
  /*
   * Gather approach:
   * For each particle i, sum over all particles j to compute rho_i.
   *
   * rho_i = sum_j m_j W(|r_i - r_j|, h_i)
   *
   * Since each OpenMP thread only updates particles[i],
   * this loop is safe to parallelize.
   */

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {

    double local_rho = 0.0;
    double drhodh = 0.0;

    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {

      Particle *p_j = &sph->particles[j];

      // Calculate 3D distance
      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;
      double dz = p_i->z - p_j->z;

      double r = sqrt(dx * dx + dy * dy + dz * dz);

      /*
       * Gather approach:
       * use particle i's smoothing length h_i.
       */
      if (r <= p_i->h) {

        double W = 0.0;
        double dWdr = 0.0;
        double dWdh = 0.0;

        cubic_spline_kernel_3d(r, p_i->h, &W, &dWdr, &dWdh);

        // rho_i = sum_j m_j W_ij(h_i)
        local_rho += p_j->mass * W;

        // drho_i / dh_i = sum_j m_j dW_ij / dh_i
        drhodh += p_j->mass * dWdh;
      }
    }

    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}

void compute_density_xreflective_yperiodic(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
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
        cubic_spline_kernel(r, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      // Left mirror (reflected across x=0)
      double dx_L = p_i->x + p_j->x;
      double r_L = sqrt(dx_L * dx_L + dy * dy);
      if (r_L <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel(r_L, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      // Right mirror (reflected across x=box_size_x)
      double dx_R = p_i->x + p_j->x - 2.0 * sph->box_size_x;
      double r_R = sqrt(dx_R * dx_R + dy * dy);
      if (r_R <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel(r_R, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }
    }

    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}

void compute_density_xperiodic_yperiodic(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {
    double local_rho = 0.0;
    double drhodh = 0.0;
    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j];

      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;

      // X-Periodic
      if (dx > 0.5 * sph->box_size_x)
        dx -= sph->box_size_x;
      else if (dx < -0.5 * sph->box_size_x)
        dx += sph->box_size_x;

      // Y-Periodic
      if (dy > 0.5 * sph->box_size_y)
        dy -= sph->box_size_y;
      else if (dy < -0.5 * sph->box_size_y)
        dy += sph->box_size_y;

      double r = sqrt(dx * dx + dy * dy);

      if (r <= p_i->h) {
        double W = 0.0, dWdr = 0.0, dWdh = 0.0;
        cubic_spline_kernel_2d(r, p_i->h, &W, &dWdr, &dWdh);
        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }
    }

    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}

void compute_density_xreflective_yzperiodic_3d(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {

    double local_rho = 0.0;
    double drhodh = 0.0;

    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {

      Particle *p_j = &sph->particles[j];

      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;
      double dz = p_i->z - p_j->z;

      /*
       * Y-periodic minimum image
       */
      if (dy > 0.5 * sph->box_size_y) {
        dy -= sph->box_size_y;
      }
      if (dy < -0.5 * sph->box_size_y) {
        dy += sph->box_size_y;
      }

      /*
       * Z-periodic minimum image
       */
      if (dz > 0.5 * sph->box_size_z) {
        dz -= sph->box_size_z;
      }
      if (dz < -0.5 * sph->box_size_z) {
        dz += sph->box_size_z;
      }

      /*
       * Direct contribution
       */
      double r = sqrt(dx * dx + dy * dy + dz * dz);

      if (r <= p_i->h) {
        double W = 0.0;
        double dWdr = 0.0;
        double dWdh = 0.0;

        /*
         * This kernel must be the 3D cubic spline kernel.
         */
        cubic_spline_kernel(r, p_i->h, &W, &dWdr, &dWdh);

        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      /*
       * Left reflective image across x = 0
       *
       * x_j,ghost = -x_j
       * dx_L = x_i - x_j,ghost = x_i + x_j
       */
      double dx_L = p_i->x + p_j->x;
      double r_L = sqrt(dx_L * dx_L + dy * dy + dz * dz);

      if (r_L <= p_i->h) {
        double W = 0.0;
        double dWdr = 0.0;
        double dWdh = 0.0;

        cubic_spline_kernel(r_L, p_i->h, &W, &dWdr, &dWdh);

        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }

      /*
       * Right reflective image across x = box_size_x
       *
       * x_j,ghost = 2 * box_size_x - x_j
       * dx_R = x_i - x_j,ghost
       *      = x_i + x_j - 2 * box_size_x
       */
      double dx_R = p_i->x + p_j->x - 2.0 * sph->box_size_x;
      double r_R = sqrt(dx_R * dx_R + dy * dy + dz * dz);

      if (r_R <= p_i->h) {
        double W = 0.0;
        double dWdr = 0.0;
        double dWdh = 0.0;

        cubic_spline_kernel_3d(r_R, p_i->h, &W, &dWdr, &dWdh);

        local_rho += p_j->mass * W;
        drhodh += p_j->mass * dWdh;
      }
    }

    p_i->rho = local_rho;
    p_i->drho_dh = drhodh;
  }
}