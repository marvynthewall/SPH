#include "sph_system.h"

void allocate_sph_system(SPHSystem *sph, int N)
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
    sph->gamma = GAMMA;
    sph->dim = 2;

    sph->time  = 0.0;
    sph->dt    = 0.0;
    sph->t_end = 0.0;
    sph->cfl   = 0.25;
    init_sph_parameter(sph);

    sph->cell_size = 0.0;
    sph->num_cells_x = 0;
    sph->num_cells_y = 0;
    sph->total_cells = 0;
    sph->head = NULL;
    sph->next = (int *)malloc(N * sizeof(int));


    size_t bytes = N * sizeof(Particle);
    sph->particles = (Particle *)malloc(bytes);
    if (sph->particles == NULL) {
        fprintf(stderr, "Error: memory allocation failed.\n");
        sph->N = 0;
        exit(EXIT_FAILURE);
    }

#ifdef USE_GPU
    cudaError_t err = cudaMalloc((void**)&sph->d_particles, bytes);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: GPU memory allocation failed: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
#endif

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
        sph->particles[i].rho      = 0.0;
        sph->particles[i].pressure = 0.0;
        sph->particles[i].u        = 0.0;
        sph->particles[i].h        = 0.0;

        sph->particles[i].cs = 0.0;
    }
}

#ifdef USE_GPU
void copy_particles_H2D(SPHSystem *sph) {
    size_t bytes = sph->N * sizeof(Particle);
    cudaMemcpy(sph->d_particles, sph->particles, bytes, cudaMemcpyHostToDevice);
}

void copy_particles_D2H(SPHSystem *sph) {
    size_t bytes = sph->N * sizeof(Particle);
    cudaMemcpy(sph->particles, sph->d_particles, bytes, cudaMemcpyDeviceToHost);
}
#endif

void free_sph_system(SPHSystem *sph)
{
    if (sph == NULL) {
        return;
    }

    if (sph->particles){
        free(sph->particles);
        sph->particles = NULL;
    }
#ifdef USE_GPU
    if (sph->d_particles) { 
        cudaFree(sph->d_particles); 
        sph->d_particles = NULL; 
    }
#endif
    sph->N = 0;

    sph->time  = 0.0;
    sph->dt    = 0.0;
    sph->t_end = 0.0;
    sph->cfl   = 0.0;
    sph->epsilon = 0.0;
    sph->alpha = 0.0;
    sph->beta = 0.0;

    free(sph->head); 
    sph->head = NULL;
    free(sph->next); 
    sph->next = NULL;
}
