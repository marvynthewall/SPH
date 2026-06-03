#include "density.cuh"

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
    /*Newton iteration: Repeatedly update h until the maximum relative change of h is small enough.*/
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];
        /* 2D adaptive h target:
         *     rho_i * h_i^2 = eta^2 * m_i
         * Define:
         *     C = eta^2 * m_i*/
        double C = eta * eta * p_i->mass;

        double h = p_i->h;
        double rho = p_i->rho;
        double drhodh = p_i->drho_dh;
         /* Newton method solves:
         *     F(h) = rho(h) * h^2 - C = 0
         * Its derivative is:
         *     dF/dh = (drho/dh) * h^2 + 2 * rho * h*/
        double F = rho * h * h - C;
        double dFdh = drhodh * h * h + 2.0 * rho * h;

        /* Avoid division by a very small number.*/
        if (fabs(dFdh) < 1.0e-14) {continue;}

        double h_new = h - F / dFdh;

        /*Reject invalid smoothing length.*/
        if (!isfinite(h_new) || h_new <= 0.0) {continue;}

        /*Compute relative change of h for convergence check.*/
        double change = fabs(h_new - h) / h;

        if (change > max_change) {max_change = change;}

        p_i->h = h_new;
    }
    /*
     * If the largest relative change is smaller than tol,
     * the adaptive h iteration is considered converged.
     */
    if (max_change < tol) {
      // printf("adaptive h 2D converged at iter: %d, eta = %.3f\n", iter, eta);
      break;
    }
  }
}
