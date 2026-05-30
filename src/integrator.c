#include "integrator.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double compute_timestep(SPHSystem2D *sph)
{
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;

    for (int i = 0; i < sph->N; i++) {
        double cs = sph->particles[i].cs;
        double h  = sph->particles[i].h;

        if (h > 0.0 && cs > 0.0) {
            double dt_i = cfl * h / cs;
            dt_min = min(dt_min, dt_i);
        }
    }

    sph->dt = dt_min;
    return dt_min;
}


double compute_timestep_signal_velocity(SPHSystem2D *sph)
{
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;

    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        double h_i = p_i->h;
        double vmax_i = p_i->cs;

        for (int j = 0; j < sph->N; j++) {

            if (i == j) continue;

            Particle *p_j = &sph->particles[j];

            double dx = p_i->x - p_j->x;
            double dy = p_i->y - p_j->y;

            double r = sqrt(dx * dx + dy * dy);

            if (r < 1.0e-12) continue;

            /*
             * If your kernel uses q = r / (2h),
             * the support radius is 2h.
             */
            if (r > 2.0 * h_i) continue;

            double dvx = p_i->vx - p_j->vx;
            double dvy = p_i->vy - p_j->vy;

            double vij_dot_rij = dvx * dx + dvy * dy;

            /*
             * wij = (v_i - v_j) dot (r_i - r_j) / |r_i - r_j|
             * wij < 0 means particles are approaching.
             */
            double wij = vij_dot_rij / r;

            /*
             * Signal velocity.
             * Basic part: c_i + c_j
             * Compression correction is added only when wij < 0.
             */
            double vsig_ij = p_i->cs + p_j->cs;

            if (wij < 0.0) {
                vsig_ij -= 3.0 * wij;
            }

            if (vsig_ij > vmax_i) {
                vmax_i = vsig_ij;
            }
        }

        if (h_i > 0.0 && vmax_i > 0.0) {
            double dt_i = cfl * h_i / vmax_i;

            if (dt_i < dt_min) {
                dt_min = dt_i;
            }
        }
    }

    if (dt_min == DBL_MAX || !isfinite(dt_min) || dt_min <= 0.0) {
        fprintf(stderr,
                "Error: invalid timestep in compute_timestep_signal_velocity. dt=%e\n",
                dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = dt_min;
    return dt_min;
}


double step_euler(
    SPHSystem2D *sph,
    double (*calculate_timep_step)(SPHSystem2D *),
    void (*compute_forces)(SPHSystem2D *)
)
{
    double dt = calculate_timep_step(sph);

    for (int i = 0; i < sph->N; i++) {

        // use old velocity to update position
        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;

        // then update velocity using old acceleration
        sph->particles[i].vx += sph->particles[i].ax * dt;
        sph->particles[i].vy += sph->particles[i].ay * dt;

        // update internal energy using old dudt
        sph->particles[i].u += sph->particles[i].dudt * dt;

        if (sph->particles[i].u < 1e-10) {
            sph->particles[i].u = 1e-10;
        }
    }

    // update hydrodynamic quantities for next step
    compute_density(sph);
    compute_pressure_soundspeed_factor(sph);
    compute_forces(sph);

    sph->time += dt;

    return dt;
}

double step_leapfrog_kdk(
    SPHSystem2D *sph,
    double (*calculate_time_step)(SPHSystem2D *),
    void (*compute_forces)(SPHSystem2D *)
)
{
    double dt = calculate_time_step(sph);

    // Kick: half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

        sph->particles[i].u  += 0.5 * sph->particles[i].dudt * dt;

        if (sph->particles[i].u < 1e-10) {
            sph->particles[i].u = 1e-10;
        }
    }

    // Drift: full-step position update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;
    }

    sph->time += dt;

    // Update hydrodynamic quantities at new position
    compute_density(sph);
    compute_pressure_soundspeed_factor(sph);
    compute_forces(sph);

    // Kick: another half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

        sph->particles[i].u  += 0.5 * sph->particles[i].dudt * dt;

        if (sph->particles[i].u < 1e-10) {
            sph->particles[i].u = 1e-10;
        }
    }
    return dt;
}