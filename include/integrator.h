#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sph_system.h"
# define min(a, b) (((a) < (b)) ? (a) : (b))


double compute_timestep(SPHSystem2D * sph);
double compute_timestep_signal_velocity(SPHSystem2D *sph);

double step_euler(SPHSystem2D *sph, double (*calculate_timep_step)(SPHSystem2D *) , void (*compute_forces)(SPHSystem2D *));
double step_leapfrog_kdk(SPHSystem2D *sph, double (*calculate_timep_step)(SPHSystem2D *) , void (*compute_forces)(SPHSystem2D *));

double step_euler_xreflective_yperiodic(SPHSystem2D *sph,
                                        double (*calculate_time_step)(SPHSystem2D *),
                                        void (*compute_forces)(SPHSystem2D *));

double step_leapfrog_kdk_xreflective_yperiodic(SPHSystem2D *sph,
                                               double (*calculate_time_step)(SPHSystem2D *),
                                               void (*compute_forces)(SPHSystem2D *));

double step_leapfrog_kdk_xperiodic_yperiodic(SPHSystem2D *sph,
                                             double (*calculate_time_step)(SPHSystem2D *),
                                             void (*compute_forces)(SPHSystem2D *));

#endif