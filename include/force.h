#ifndef FORCE_H
#define FORCE_H

#include "sph_system.h"

void compute_pressure_soundspeed_factor(SPHSystem *sph);
void compute_force(SPHSystem *sph);
void compute_force_1d_xreflective(SPHSystem *sph);
void compute_force_xreflective_yperiodic(SPHSystem *sph);
void compute_force_xperiodic_yperiodic(SPHSystem *sph);
void compute_force_xreflective_yperiodic_celllist(SPHSystem *sph);

void compute_force_3d(SPHSystem *sph);
void compute_force_xreflective_yperiodic_zperiodic_3d(SPHSystem *sph);

#endif
