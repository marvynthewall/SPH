#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* Particle Structures - Including 2D SPH fields */
typedef struct {

    int id;

    // position
    double x, y;

    // velocity
    double vx, vy;

    // acceleration
    double ax, ay;

    // SPH particle properties
    double mass;       // m
    double density;    // rho
    double pressure;   // P
    double u;          // specific internal energy
    double h;          // smoothing length

    // useful derived quantities
    double cs;         // sound speed

} Particle;


/* SPH System Structure */
typedef struct {

    int N;

    Particle *particles;

    // time control
    double time;
    double dt;
    double t_end;
    double cfl;

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

#endif