#include "../include/integrator.h"

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

double step_euler(SPHSystem2D *sph)
{
    double dt = compute_timestep(sph);

    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].vx += sph->particles[i].ax * dt;
        sph->particles[i].vy += sph->particles[i].ay * dt;

        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;
    }
    sph->time += dt;
    return sph->time;
}


double step_leapfrog_kdk(SPHSystem2D *sph, void (*compute_forces)(SPHSystem2D *))
{
    double dt = compute_timestep(sph);

    // Kick: half-step velocity update
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;
    }

    // Drift: full-step position update
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;
    }

    sph->time += dt;

    compute_forces(sph);
    // Kick: another half-step velocity update
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;
    }


    return dt;
}