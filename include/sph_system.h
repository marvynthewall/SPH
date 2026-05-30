#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

// Physical constants
// can be put into constants.h file
// for most case in the Astrophysics, or shockwave in the ISM (Interstellar Medium)
// only consider monatomic gas gamma = 1.667
// if for diatomic gas gamma = 1.4

// for 2D, ideal monatonic gas gamma = (f+2) / f = 4 / 2 = 2.0
#define GAMMA 1.666666666666667


/* Particle Structures - Including 2D SPH fields */
typedef struct {

    int id;

    // position
    double x, y;

    // velocity
    double vx, vy;

    // acceleration
    double ax, ay;

    // specific internal energy evolution
    double dudt;

    // SPH particle properties
    double mass;       // m
    double rho;        // rho
    double pressure;   // P
    double u;          // specific internal energy
    double h;          // smoothing length

    // useful derived quantities
    double cs;         // sound speed
    double drho_dh;     // for pressure calculation
    double factor;      // for pressure

} Particle;


/* SPH System Structure */
typedef struct {

    int N;

    Particle *particles;

    // time control
    double time;
    double dt;
    double t_end;
    // 0.1 - 0.3 from paper
    double cfl;

    // for Viscositvy 
    double epsilon;
    double alpha;
    double beta;
} SPHSystem2D;


/* memory management */
void allocate_sph_system(SPHSystem2D *sph, int N);
void free_sph_system(SPHSystem2D *sph);


/* module headers */
#include "density.h"
#include "force.h"
#include "init.h"
#include "integrator.h"
#include "io.h"
#include "kernel.h"
#include "constants.h"

#endif
