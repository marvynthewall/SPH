#pragma once
#include "sph_system.h"

void build_cell_list(SPHSystem *sph);

double compute_timestep(SPHSystem *sph);
double compute_timestep_signal_velocity(SPHSystem *sph);
double compute_timestep_signal_velocity_3d(SPHSystem *sph);

double step_euler(SPHSystem *sph, double (*calculate_timep_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *));
double step_leapfrog_kdk(SPHSystem *sph,
        double (*calculate_timep_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *));

double
step_euler_xreflective_yperiodic(SPHSystem *sph,
                                 double (*calculate_time_step)(SPHSystem *),
                                 void (*compute_forces)(SPHSystem *));
double step_euler_xreflective_yzperiodic_3d(
    SPHSystem *sph, double (*calculate_timep_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *));

double step_leapfrog_kdk_xreflective_yperiodic(
    SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *));
double step_leapfrog_kdk_xperiodic_yperiodic(
    SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *));
double step_leapfrog_kdk_xreflective_yzperiodic_3d(
    SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *));
double
step_leapfrog_kdk_1d_xreflective(SPHSystem *sph,
                                 double (*calculate_time_step)(SPHSystem *),
                                 void (*compute_forces)(SPHSystem *));


// GPU
double step_leapfrog_kdk_xreflective_yperiodic_gpu(
        SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
        void (*compute_forces)(SPHSystem *));

double compute_timestep_signal_velocity_gpu(SPHSystem *sph);
