#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <errno.h>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

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
    double x, y, z;

    // velocity
    double vx, vy, vz;

    // acceleration
    double ax, ay, az;

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
    int dim;
    double gamma;

    int N;

    Particle *particles;

#ifdef __CUDACC__
    // This is the copy in GPU
    // Actually all the computation run this
    // So the CPU version is actually the copy of this GPU
    Particle *d_particles;
#endif

    // time control
    double time;
    double dt;
    double t_end;

    // domain size
    double box_size_x;
    double box_size_y;
    double box_size_z;

    // 0.1 - 0.3 from paper
    double cfl;

    // for Viscositvy 
    double epsilon;
    double alpha;
    double beta;

    // --- Cell-Linked List variables ---
    double cell_size;     // size (>= max_h)
    int num_cells_x;      // X cells
    int num_cells_y;      // Y cells
    int num_cells_z;      // Z cells
    int total_cells;      // = num_cells_x * num_cells_y

    int *head;            // array size = total_cells
    int *next;            // array size = N
#ifdef __CUDACC__
    int *d_head;
    int *d_next;
#endif
} SPHSystem;


/* memory management */
void allocate_sph_system(SPHSystem *sph, int N);
void free_sph_system(SPHSystem *sph);

#ifdef __CUDACC__
void copy_particles_H2D(SPHSystem *sph);
void copy_particles_D2H(SPHSystem *sph);
#endif

#endif
