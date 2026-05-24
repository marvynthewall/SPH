#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// temp
typedef struct {
    int N;

    double *x;
    double *y;

    double *vx;
    double *vy;

    double *ax;
    double *ay;

    double *m;
    double *rho;
    double *P;
    double *u;
    double *h;

} SPHSystem2D;

void allocate_sph_system(SPHSystem2D *sph, int N);
void free_sph_system(SPHSystem2D *sph);


# include "density.h"
# include "force.h"
# include "init.h"
# include "integrator.h"
# include "io.h"
# include "kernel.h"


#endif