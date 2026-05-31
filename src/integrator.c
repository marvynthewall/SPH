#include "integrator.h"

#ifdef _OPENMP
#include <omp.h>
#endif


double compute_timestep(SPHSystem *sph)
{
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:dt_min)  schedule(dynamic)
#endif
    for (int i = 0; i < sph->N; i++) {
        double cs = sph->particles[i].cs;
        double h  = sph->particles[i].h;

        if (h > 0.0 && cs > 0.0) {
            double dt_i = cfl * h / cs;
            dt_min = min(dt_min, dt_i);
        }
    }

    if (dt_min == DBL_MAX || !isfinite(dt_min) || dt_min <= 0.0) {
        fprintf(stderr,
                "Error: invalid timestep in compute_timestep. dt=%e\n",
                dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = dt_min;
    return dt_min;
}


double compute_timestep_signal_velocity(SPHSystem *sph)
{
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;

#ifdef _OPENMP
#pragma omp parallel for reduction(min:dt_min) schedule(dynamic)
#endif
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

            /* If your kernel uses q = r / (h), the support radius is h. */
            if (r > h_i) continue;

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
            dt_min = min(dt_min, dt_i);
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
    SPHSystem *sph,
    double (*calculate_timep_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *)
)
{
    double dt = calculate_timep_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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
    SPHSystem *sph,
    double (*calculate_time_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *)
)
{
    double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Kick: half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

        sph->particles[i].u  += 0.5 * sph->particles[i].dudt * dt;

        if (sph->particles[i].u < 1e-10) {
            sph->particles[i].u = 1e-10;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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



double step_euler_xreflective_yperiodic(SPHSystem *sph,
                  double (*calculate_timep_step)(SPHSystem *),
                  void (*compute_forces)(SPHSystem *)) {
  double dt = calculate_timep_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    // use old velocity to update position
    sph->particles[i].x += sph->particles[i].vx * dt;
    sph->particles[i].y += sph->particles[i].vy * dt;

    // Y-Periodic
    if (sph->particles[i].y >= sph->box_size_y)
      sph->particles[i].y -= sph->box_size_y;
    if (sph->particles[i].y < 0.0)
      sph->particles[i].y += sph->box_size_y;

    // X-Reflective
    if (sph->particles[i].x < 0.0) {
      sph->particles[i].x = -sph->particles[i].x;
      sph->particles[i].vx = -sph->particles[i].vx;
    }
    if (sph->particles[i].x > sph->box_size_x) {
      sph->particles[i].x = 2.0 * sph->box_size_x - sph->particles[i].x;
      sph->particles[i].vx = -sph->particles[i].vx;
    }

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
  update_adaptive_h_2d(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic);
  compute_density_xreflective_yperiodic(sph);
  compute_pressure_soundspeed_factor(sph);
  compute_forces(sph);

  sph->time += dt;

  return dt;
}

double step_leapfrog_kdk_xreflective_yperiodic(SPHSystem *sph,
                         double (*calculate_time_step)(SPHSystem *),
                         void (*compute_forces)(SPHSystem *)) {
  double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Kick: half-step velocity and internal energy update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;

    if (sph->particles[i].u < 1e-10) {
      sph->particles[i].u = 1e-10;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Drift: full-step position update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].x += sph->particles[i].vx * dt;
    sph->particles[i].y += sph->particles[i].vy * dt;

    // Y-Periodic
    if (sph->particles[i].y >= sph->box_size_y)
      sph->particles[i].y -= sph->box_size_y;
    if (sph->particles[i].y < 0.0)
      sph->particles[i].y += sph->box_size_y;

    // X-Reflective
    if (sph->particles[i].x < 0.0) {
      sph->particles[i].x = -sph->particles[i].x;
      sph->particles[i].vx = -sph->particles[i].vx;
    }
    if (sph->particles[i].x > sph->box_size_x) {
      sph->particles[i].x = 2.0 * sph->box_size_x - sph->particles[i].x;
      sph->particles[i].vx = -sph->particles[i].vx;
    }
  }

  sph->time += dt;

  // Update hydrodynamic quantities at new position
  update_adaptive_h_2d(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic);
  compute_density_xreflective_yperiodic(sph);
  compute_pressure_soundspeed_factor(sph);
  compute_forces(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Kick: another half-step velocity and internal energy update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;

    if (sph->particles[i].u < 1e-10) {
      sph->particles[i].u = 1e-10;
    }
  }
  return dt;
}

double step_leapfrog_kdk_xperiodic_yperiodic(SPHSystem *sph,
                         double (*calculate_time_step)(SPHSystem *),
                         void (*compute_forces)(SPHSystem *)) {
  double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Kick: half-step velocity and internal energy update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;

    if (sph->particles[i].u < 1e-10) {
      sph->particles[i].u = 1e-10;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Drift: full-step position update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].x += sph->particles[i].vx * dt;
    sph->particles[i].y += sph->particles[i].vy * dt;

    // X-Periodic
    if (sph->particles[i].x >= sph->box_size_x)
      sph->particles[i].x -= sph->box_size_x;
    if (sph->particles[i].x < 0.0)
      sph->particles[i].x += sph->box_size_x;

    // Y-Periodic
    if (sph->particles[i].y >= sph->box_size_y)
      sph->particles[i].y -= sph->box_size_y;
    if (sph->particles[i].y < 0.0)
      sph->particles[i].y += sph->box_size_y;
  }

  sph->time += dt;

  // Update hydrodynamic quantities at new position
  update_adaptive_h_2d(sph, 20, 1e-4, 1.8, compute_density_xperiodic_yperiodic);
  compute_density_xperiodic_yperiodic(sph);
  compute_pressure_soundspeed_factor(sph);
  compute_forces(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Kick: another half-step velocity and internal energy update
  for (int i = 0; i < sph->N; i++) {

    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;

    if (sph->particles[i].u < 1e-10) {
      sph->particles[i].u = 1e-10;
    }
  }
  return dt;
}
