#include "sph_system.h"
#include "density.h"
#include "init.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void update_adaptive_h(SPHSystem2D *sph, int max_iter, double tol)
{
    const double eta = 2.3;

    if (sph == NULL || sph->particles == NULL) {
        fprintf(stderr, "Error: invalid SPH system in update_adaptive_h.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];
        double C = eta * eta * p_i->mass;

        for (int iter = 0; iter < max_iter; iter++) {
            compute_density(sph);

            double h = p_i->h;
            double rho = p_i->rho;
            double drhodh = p_i->drho_dh;

            double F = rho * h * h - C;
            double dFdh = drhodh * h * h + 2.0 * rho * h;

            if (fabs(dFdh) < 1.0e-14) {
                fprintf(stderr, "Warning: dFdh too small at particle %d\n", i);
                break;
            }

            double h_new = h - F / dFdh;

            if (!isfinite(h_new) || h_new <= 0.0) {
                fprintf(stderr,
                        "Warning: invalid h_new at particle %d: h_new=%e\n", i, h_new);
                break;
            }

            p_i->h = h_new;

            if (fabs(h_new - h) / h < tol) {
                break;
            }
        }
    }

    /*
     * Important:
     * Recompute density once using the final updated h values.
     */
    compute_density(sph);
}

void check_adaptive_h(SPHSystem2D *sph, double eta, double tol)
{
    if (sph == NULL || sph->particles == NULL) {
        fprintf(stderr, "Error: invalid SPH system in check_adaptive_h.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        double C = eta * eta * p_i->mass;
        double residual =
            fabs(p_i->rho * p_i->h * p_i->h - C) / C;

        printf("i=%d h=%e rho=%e rho*h^2=%e C=%e residual=%e\n",
               i,
               p_i->h,
               p_i->rho,
               p_i->rho * p_i->h * p_i->h,
               C,
               residual);

        if (residual > tol) {
            fprintf(stderr, "Warning: particle %d not converged, residual=%e\n", i, residual);
        }
    }
}


int main(int argc, char *argv[]) {
    SPHSystem2D sph;

    int nx = 20; 
    int ny = 20;
    int N  = nx * ny;

    allocate_sph_system(&sph, N);
    init_uniform_box(&sph, nx, ny);

    update_adaptive_h(&sph, 100, 1.0e-4);
    check_adaptive_h(&sph, 2.3, 1.0e-4);
    free_sph_system(&sph);
    
    return 0;
}