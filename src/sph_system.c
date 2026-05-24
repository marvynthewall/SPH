#include "sph_system.h"

void allocate_sph_system(SPHSystem2D *sph, int N)
{
    if (sph == NULL) {
        fprintf(stderr, "Error: sph is NULL.\n");
        exit(EXIT_FAILURE);
    }

    if (N <= 0) {
        fprintf(stderr,
                "Error: N must be positive. N = %d\n",
                N);
        exit(EXIT_FAILURE);
    }

    sph->N = N;

    sph->time  = 0.0;
    sph->dt    = 0.0;
    sph->t_end = 0.0;
    sph->cfl   = 0.25;

    sph->particles = (Particle *)malloc(N * sizeof(Particle));

    if (sph->particles == NULL) {
        fprintf(stderr, "Error: memory allocation failed.\n");
        sph->N = 0;
        exit(EXIT_FAILURE);
    }

    // initialize particles
    for (int i = 0; i < N; i++) {

        sph->particles[i].id = i;

        // position
        sph->particles[i].x = 0.0;
        sph->particles[i].y = 0.0;

        // velocity
        sph->particles[i].vx = 0.0;
        sph->particles[i].vy = 0.0;

        // acceleration
        sph->particles[i].ax = 0.0;
        sph->particles[i].ay = 0.0;

        // SPH quantities
        sph->particles[i].mass     = 0.0;
        sph->particles[i].density  = 0.0;
        sph->particles[i].pressure = 0.0;
        sph->particles[i].u        = 0.0;
        sph->particles[i].h        = 0.0;

        sph->particles[i].cs = 0.0;
    }
}

void free_sph_system(SPHSystem2D *sph)
{
    if (sph == NULL) {
        return;
    }

    free(sph->particles);

    sph->particles = NULL;

    sph->particles = NULL;
    sph->N = 0;

    sph->time  = 0.0;
    sph->dt    = 0.0;
    sph->t_end = 0.0;
    sph->cfl   = 0.0;
}