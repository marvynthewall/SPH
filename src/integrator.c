#include "sph_all.h"

double compute_timestep(SPHSystem *sph) {
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;
#ifdef _OPENMP
#pragma omp parallel for reduction(min : dt_min) schedule(dynamic)
#endif
    for (int i = 0; i < sph->N; i++) {
        double cs = sph->particles[i].cs;
        double h = sph->particles[i].h;

        if (h > 0.0 && cs > 0.0) {
            double dt_i = cfl * h / cs;
            dt_min = fmin(dt_min, dt_i);
        }
    }

    if (dt_min == DBL_MAX || !isfinite(dt_min) || dt_min <= 0.0) {
        fprintf(stderr, "Error: invalid timestep in compute_timestep. dt=%e\n",
                dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = dt_min;
    return dt_min;
}

double compute_timestep_signal_velocity(SPHSystem *sph) {
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;

#ifdef _OPENMP
#pragma omp parallel for reduction(min : dt_min) schedule(dynamic)
#endif
    for (int i = 0; i < sph->N; i++) {
        Particle *p_i = &sph->particles[i];

        double h_i = p_i->h;
        double vmax_i = p_i->cs;

        for (int j = 0; j < sph->N; j++) {

            if (i == j)
                continue;

            Particle *p_j = &sph->particles[j];

            double dx = p_i->x - p_j->x;
            double dy = p_i->y - p_j->y;

            double r = sqrt(dx * dx + dy * dy);

            if (r < 1.0e-12)
                continue;

            /* If your kernel uses q = r / (h), the support radius is h. */
            if (r > h_i)
                continue;

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

            vsig_ij -= (wij < 0.0) ? 3.0 * wij : 0.0;

            vmax_i = (vsig_ij > vmax_i) ? vsig_ij : vmax_i;
        }

        if (h_i > 0.0 && vmax_i > 0.0) {
            double dt_i = cfl * h_i / vmax_i;
            if(dt_i > 0.01)dt_i = 0.01;
            dt_min = (dt_i < dt_min) ? dt_i : dt_min;
        }
    }

    if (dt_min == DBL_MAX || !isfinite(dt_min) || dt_min <= 0.0) {
        fprintf(
                stderr,
                "Error: invalid timestep in compute_timestep_signal_velocity. dt=%e\n",
                dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = dt_min;
    return dt_min;
}

double compute_timestep_signal_velocity_3d(SPHSystem *sph)
{
    double dt_min = DBL_MAX;
    double cfl = sph->cfl;

#ifdef _OPENMP
#pragma omp parallel for reduction(min : dt_min) schedule(static)
#endif
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        double h_i = p_i->h;
        double vmax_i = p_i->cs;

        for (int j = 0; j < sph->N; j++) {

            if (i == j)
                continue;

            Particle *p_j = &sph->particles[j];

            double dx = p_i->x - p_j->x;

            double dy = p_i->y - p_j->y;
            // Cache the box size in a local variable (CPU register) to optimize memory access
            double box_y = sph->box_size_y;
            // Minimum image convention for periodic distance via branchless ternary operators
            dy -= (dy >  0.5 * box_y) ? box_y : 0.0;
            dy += (dy < -0.5 * box_y) ? box_y : 0.0;

            double dz = p_i->z - p_j->z;
            // Cache the box size in a local variable (CPU register) to optimize memory access
            double box_z = sph->box_size_z;
            // Minimum image convention for periodic distance via branchless ternary operators
            dz -= (dz >  0.5 * box_z) ? box_z : 0.0;
            dz += (dz < -0.5 * box_z) ? box_z : 0.0;

            double r = sqrt(dx * dx + dy * dy + dz * dz);

            if (r < 1.0e-12)
                continue;

            if (r > h_i)
                continue;

            double dvx = p_i->vx - p_j->vx;
            double dvy = p_i->vy - p_j->vy;
            double dvz = p_i->vz - p_j->vz;

            double vij_dot_rij = dvx * dx + dvy * dy + dvz * dz;
            double wij = vij_dot_rij / r;

            double vsig_ij = p_i->cs + p_j->cs;

            vsig_ij -= (wij < 0.0) ? 3.0 * wij : 0.0;

            vmax_i = (vsig_ij > vmax_i) ? vsig_ij : vmax_i;
        }

        if (h_i > 0.0 && vmax_i > 0.0) {
            double dt_i = cfl * h_i / vmax_i;
            // 4. Update global minimum timestep (Hardware-accelerated fmin)
            dt_min = (dt_i < dt_min) ? dt_i : dt_min;
        }
    }

    if (dt_min == DBL_MAX || !isfinite(dt_min) || dt_min <= 0.0) {
        fprintf(stderr,
                "Error: invalid timestep in compute_timestep_signal_velocity_3d. dt=%e\n",
                dt_min);
        exit(EXIT_FAILURE);
    }

    sph->dt = dt_min;
    return dt_min;
}

double step_euler(SPHSystem *sph, double (*calculate_timep_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *)) {
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

double step_leapfrog_kdk(SPHSystem *sph,
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

        sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;

        if (sph->particles[i].u < 1e-10) {
            sph->particles[i].u = 1e-10;
        }
    }
    return dt;
}

double
step_euler_xreflective_yperiodic(SPHSystem *sph,
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
    update_adaptive_h(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic);
    compute_density_xreflective_yperiodic(sph);
    compute_pressure_soundspeed_factor(sph);
    compute_forces(sph);

    sph->time += dt;

    return dt;
}

double
step_euler_xreflective_yzperiodic_3d(SPHSystem *sph,
        double (*calculate_time_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *)) {
    double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < sph->N; i++) {

        Particle *p = &sph->particles[i];

        /*
         * Explicit Euler position update:
         * x^{n+1} = x^n + v^n dt
         */
        p->x += p->vx * dt;
        p->y += p->vy * dt;
        p->z += p->vz * dt;

        /*
         * Explicit Euler velocity update:
         * v^{n+1} = v^n + a^n dt
         */
        p->vx += p->ax * dt;
        p->vy += p->ay * dt;
        p->vz += p->az * dt;

        /*
         * Internal energy update:
         * u^{n+1} = u^n + dudt^n dt
         */
        p->u += p->dudt * dt;

        if (p->u < 1.0e-10) {
            p->u = 1.0e-10;
        }

        /* Y-periodic boundary */
        while (p->y >= sph->box_size_y) {
            p->y -= sph->box_size_y;
        }
        while (p->y < 0.0) {
            p->y += sph->box_size_y;
        }

        /* Z-periodic boundary */
        while (p->z >= sph->box_size_z) {
            p->z -= sph->box_size_z;
        }
        while (p->z < 0.0) {
            p->z += sph->box_size_z;
        }

        /* X-reflective boundary */
        if (p->x < 0.0) {
            p->x = -p->x;
            p->vx = -p->vx;
        }
        if (p->x > sph->box_size_x) {
            p->x = 2.0 * sph->box_size_x - p->x;
            p->vx = -p->vx;
        }
    }

    /*
     * Update hydrodynamic quantities for next step.
     * These functions must also be 3D and boundary-consistent.
     */
    update_adaptive_h_3d(sph, 20, 1.0e-4, 2.3,
            compute_density_xreflective_yzperiodic_3d);

    compute_density_xreflective_yzperiodic_3d(sph);
    compute_pressure_soundspeed_factor(sph);
    compute_forces(sph);

    sph->time += dt;
    sph->dt = dt;

    return dt;
}

double step_leapfrog_kdk_xreflective_yperiodic(
        SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *)) {
    double dt = calculate_time_step(sph);
#ifdef _OPENMP
    double t1 = omp_get_wtime();
#else
    struct timeval t1, t2, t3, t4, t5, t6, t7;
    gettimeofday(&t1, NULL);
#endif


#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Kick: half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

        sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;
        sph->particles[i].u = fmax(1e-10, sph->particles[i].u);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Drift: full-step position update
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;

        // 1. Load data into local variables (CPU registers) to minimize memory access overhead
        double x     = sph->particles[i].x;
        double y     = sph->particles[i].y;
        double vx    = sph->particles[i].vx;
        double box_x = sph->box_size_x;
        double box_y = sph->box_size_y;

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
        sph->particles[i].x  = x;
        sph->particles[i].y  = y;
        sph->particles[i].vx = vx;
    }

    sph->time += dt;

#ifdef _OPENMP
    double t2 = omp_get_wtime();
#else
    gettimeofday(&t2, NULL);
#endif


    // Update hydrodynamic quantities at new position
    // update_adaptive_h(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic);
    update_adaptive_h(sph, 20, 1e-4, 2.3, compute_density_xreflective_yperiodic_celllist);



#ifdef _OPENMP
    double t3 = omp_get_wtime();
#else
    gettimeofday(&t3, NULL);
#endif
    // compute_density_xreflective_yperiodic(sph);
    compute_density_xreflective_yperiodic_celllist(sph);
#ifdef _OPENMP
    double t4 = omp_get_wtime();
#else
    gettimeofday(&t4, NULL);
#endif
    compute_pressure_soundspeed_factor(sph);
#ifdef _OPENMP
    double t5 = omp_get_wtime();
#else
    gettimeofday(&t5, NULL);
#endif
    compute_forces(sph);
#ifdef _OPENMP
    double t6 = omp_get_wtime();
#else
    gettimeofday(&t6, NULL);
#endif



#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Kick: another half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
        sph->particles[i].vy += 0.5 * sph->particles[i].ay * dt;

        sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;
        sph->particles[i].u = fmax(1e-10, sph->particles[i].u);
    }
#ifdef _OPENMP
    double t7 = omp_get_wtime();
#else
    gettimeofday(&t7, NULL);
#endif

#ifdef _OPENMP
    printf("Kick/Drift: %f s | Adapt H: %f s | Density: %f s | Pressure: %f s | Force: %f s | Kick: %f s\n", 
            t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6);
#else
    double dt21 = (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec) / 1000000.0;
    double dt32 = (t3.tv_sec-t2.tv_sec) + (t3.tv_usec-t2.tv_usec) / 1000000.0;
    double dt43 = (t4.tv_sec-t3.tv_sec) + (t4.tv_usec-t3.tv_usec) / 1000000.0;
    double dt54 = (t5.tv_sec-t4.tv_sec) + (t5.tv_usec-t4.tv_usec) / 1000000.0;
    double dt65 = (t6.tv_sec-t5.tv_sec) + (t6.tv_usec-t5.tv_usec) / 1000000.0;
    double dt76 = (t7.tv_sec-t6.tv_sec) + (t7.tv_usec-t6.tv_usec) / 1000000.0;
    printf("Kick/Drift: %f s | Adapt H: %f s | Density: %f s | Pressure: %f s | Force: %f s | Kick: %f s\n", 
            dt21, dt32, dt43, dt54, dt65, dt76);
#endif
    return dt;
}

double step_leapfrog_kdk_xperiodic_yperiodic(
        SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
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
        sph->particles[i].u = fmax(sph->particles[i].u, 1e-10);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Drift: full-step position update
    for (int i = 0; i < sph->N; i++) {

        sph->particles[i].x += sph->particles[i].vx * dt;
        sph->particles[i].y += sph->particles[i].vy * dt;

        // to register
        double x = sph->particles[i].x;
        double y = sph->particles[i].y;
        double box_x = sph->box_size_x;
        double box_y = sph->box_size_y;

        // X-Periodic
        x -= (x >= box_x) ? box_x : 0.0;
        x += (x < 0.0)    ? box_x : 0.0;

        // Y-Periodic
        y -= (y >= box_y) ? box_y : 0.0;
        y += (y < 0.0)    ? box_y : 0.0;

        // 2. go back to RAM
        sph->particles[i].x = x;
        sph->particles[i].y = y;
    }

    sph->time += dt;

    // Update hydrodynamic quantities at new position
    update_adaptive_h(sph, 20, 1e-4, 1.8, compute_density_xperiodic_yperiodic_celllist);
    compute_density_xperiodic_yperiodic_celllist(sph);
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
        sph->particles[i].u = fmax(sph->particles[i].u, 1e-10);
    }
    return dt;
}

double step_leapfrog_kdk_xreflective_yzperiodic_3d(
        SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *)) {
    double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Kick: half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        Particle *p = &sph->particles[i];

        p->vx += 0.5 * p->ax * dt;
        p->vy += 0.5 * p->ay * dt;
        p->vz += 0.5 * p->az * dt;

        p->u += 0.5 * p->dudt * dt;
        p->u = fmax(1.0e-10, p->u);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Drift: full-step position update
    for (int i = 0; i < sph->N; i++) {

        Particle *p = &sph->particles[i];

        p->x += p->vx * dt;
        p->y += p->vy * dt;
        p->z += p->vz * dt;
        // 1. Load data into local variables (CPU registers) to minimize memory access overhead
        double x     = p->x;
        double y     = p->y;
        double z     = p->z;
        double vx    = p->vx;
        double box_x = sph->box_size_x;
        double box_y = sph->box_size_y;
        double box_z = sph->box_size_z;

        // 2. Y-Periodic boundary condition (Replaced 'while' with branchless ternary operators)
        y -= (y >= box_y) ? box_y : 0.0;
        y += (y < 0.0)    ? box_y : 0.0;

        // 3. Z-Periodic boundary condition (Replaced 'while' with branchless ternary operators)
        z -= (z >= box_z) ? box_z : 0.0;
        z += (z < 0.0)    ? box_z : 0.0;

        // 4. X-Reflective boundary condition (Branchless via condition flags)
        // Handle left boundary reflection
        int hit_left = (x < 0.0);
        x  = hit_left ? -x  : x;
        vx = hit_left ? -vx : vx;

        // Handle right boundary reflection
        int hit_right = (x > box_x);
        x  = hit_right ? (2.0 * box_x - x) : x;
        vx = hit_right ? -vx : vx;

        // 5. Store updated states back to the particle structure in memory once
        p->x  = x;
        p->y  = y;
        p->z  = z;
        p->vx = vx;
    }

    sph->time += dt;

    /*
     * Update hydrodynamic quantities at the new position.
     * These functions must be 3D and boundary-consistent.
     */
    update_adaptive_h_3d(sph, 20, 1.0e-4, 2.3,compute_density_xreflective_yzperiodic_3d_celllist);
    compute_density_xreflective_yzperiodic_3d_celllist(sph);
    compute_pressure_soundspeed_factor(sph);
    compute_forces(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    // Kick: another half-step velocity and internal energy update
    for (int i = 0; i < sph->N; i++) {

        Particle *p = &sph->particles[i];

        p->vx += 0.5 * p->ax * dt;
        p->vy += 0.5 * p->ay * dt;
        p->vz += 0.5 * p->az * dt;

        p->u += 0.5 * p->dudt * dt;
        p->u = fmax(p->u, 1.0e-10);
    }

    sph->dt = dt;

    return dt;
}

double step_leapfrog_kdk_1d_xreflective(SPHSystem *sph, double (*calculate_time_step)(SPHSystem *), void (*compute_forces)(SPHSystem *)) {
  double dt = calculate_time_step(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;
    sph->particles[i].u = fmax(1e-10, sph->particles[i].u);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    double x = sph->particles[i].x + sph->particles[i].vx * dt;
    double vx = sph->particles[i].vx;
    double box_x = sph->box_size_x;

    int hit_left = (x < 0.0);
    x = hit_left ? -x : x;
    vx = hit_left ? -vx : vx;

    int hit_right = (x > box_x);
    x = hit_right ? (2.0 * box_x - x) : x;
    vx = hit_right ? -vx : vx;

    sph->particles[i].x = x;
    sph->particles[i].vx = vx;
  }

  sph->time += dt;

  update_adaptive_h(sph, 20, 1e-4, 2.3, compute_density_1d_xreflective);
  compute_density_1d_xreflective(sph);
  compute_pressure_soundspeed_factor(sph);
  compute_forces(sph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].vx += 0.5 * sph->particles[i].ax * dt;
    sph->particles[i].u += 0.5 * sph->particles[i].dudt * dt;
    sph->particles[i].u = fmax(1e-10, sph->particles[i].u);
  }

  return dt;
}
