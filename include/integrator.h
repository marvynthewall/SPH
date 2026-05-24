#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sph_system.h"
# define min(a, b) (((a) < (b)) ? (a) : (b))


double compute_timestep(SPHSystem2D * sph);


double step_euler(SPHSystem2D * sph);
double step_leapfrog_kdk(SPHSystem2D *sph, void (*compute_forces)(SPHSystem2D *));


#endif