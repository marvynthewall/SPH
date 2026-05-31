#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sph_system.h"
#define min(a, b) (((a) < (b)) ? (a) : (b))

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
double step_leapfrog_kdk_xreflective_yzperiodic_3d(
    SPHSystem *sph, double (*calculate_time_step)(SPHSystem *),
    void (*compute_forces)(SPHSystem *));

#endif
