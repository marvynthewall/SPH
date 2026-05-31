#ifndef FORCE_H
#define FORCE_H

#include "sph_system.h"

void compute_pressure_soundspeed_factor(SPHSystem *sph);
void compute_force(SPHSystem *sph);
void compute_force_xreflective_yperiodic(SPHSystem *sph);
void compute_force_3d(SPHSystem *sph);

#endif
