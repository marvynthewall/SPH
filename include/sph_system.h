#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Particle Structures - Including 2D SPH fields (Position, Velocity, Mass,
 * Density, Pressure)
 */
typedef struct {
  int id;          // particle ID
  double x, y;     // 2D position
  double vx, vy;   // 2D velocity
  double mass;     // mass (m)
  double density;  // density (rho)
  double pressure; // pressure (P)
  double h;        // smoothing length
} Particle;

#include "density.h"
#include "force.h"
#include "init.h"
#include "integrator.h"
#include "io.h"
#include "kernel.h"

#endif