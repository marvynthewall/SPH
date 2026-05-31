#ifndef FORCE_H
#define FORCE_H

#include "sph_system.h"

void compute_pressure_soundspeed_factor(SPHSystem2D *sph);
void compute_force(SPHSystem2D *sph);
void compute_force_xreflective_yperiodic(SPHSystem2D *sph);

#endif
