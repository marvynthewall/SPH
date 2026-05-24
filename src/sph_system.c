# include "sph_system.h"

// temp
void allocate_sph_system(SPHSystem2D *sph, int N)
{
    if (N <= 0) {
        fprintf(stderr, "Error: N must be positive. N = %d\n", N);
        exit(EXIT_FAILURE);
    }

    sph->N = N;

    sph->x   = (double *)malloc(N * sizeof(double));
    sph->y   = (double *)malloc(N * sizeof(double));

    sph->vx  = (double *)malloc(N * sizeof(double));
    sph->vy  = (double *)malloc(N * sizeof(double));

    sph->ax  = (double *)malloc(N * sizeof(double));
    sph->ay  = (double *)malloc(N * sizeof(double));

    sph->m   = (double *)malloc(N * sizeof(double));
    sph->rho = (double *)malloc(N * sizeof(double));
    sph->P   = (double *)malloc(N * sizeof(double));
    sph->u   = (double *)malloc(N * sizeof(double));
    sph->h   = (double *)malloc(N * sizeof(double));

    if (sph->x   == NULL || sph->y   == NULL ||
        sph->vx  == NULL || sph->vy  == NULL ||
        sph->ax  == NULL || sph->ay  == NULL ||
        sph->m   == NULL || sph->rho == NULL ||
        sph->P   == NULL || sph->u   == NULL ||
        sph->h   == NULL) {

        fprintf(stderr, "Error: memory allocation failed.\n");
        free_sph_system(sph);
        exit(EXIT_FAILURE);
    }
}


void free_sph_system(SPHSystem2D *sph)
{
    if (sph == NULL) {
        return;
    }

    free(sph->x);
    free(sph->y);

    free(sph->vx);
    free(sph->vy);

    free(sph->ax);
    free(sph->ay);

    free(sph->m);
    free(sph->rho);
    free(sph->P);
    free(sph->u);
    free(sph->h);

    sph->x   = NULL;
    sph->y   = NULL;

    sph->vx  = NULL;
    sph->vy  = NULL;

    sph->ax  = NULL;
    sph->ay  = NULL;

    sph->m   = NULL;
    sph->rho = NULL;
    sph->P   = NULL;
    sph->u   = NULL;
    sph->h   = NULL;

    sph->N = 0;
}

