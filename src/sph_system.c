#include "sph_all.h"

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

#ifdef __CUDACC__
    cudaMalloc((void**)&sph->d_particles, N * sizeof(Particle));
    cudaMalloc((void**)&sph->d_next, N * sizeof(int));
    sph->d_head = NULL; // head 會在 cell list 階段動態分配，先設為 NULL
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

#ifdef __CUDACC__
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
#ifdef __CUDACC__
    if (sph->d_particles) { 
        cudaFree(sph->d_particles); 
        sph->d_particles = NULL; 
        cudaFree(sph->d_next);
        sph->d_next = NULL;
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


void build_cell_list(SPHSystem *sph) {
    // 1. find the max h, for the safe grid size
    double max_h = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(max:max_h)
    #endif
    for (int i = 0; i < sph->N; i++) {
        if (sph->particles[i].h > max_h) {
            max_h = sph->particles[i].h;
        }
    }
    
    // prevent too many cells
    if (max_h < 1e-4) max_h = 1e-4;
    // prevent too many cells causing Integer Overflow and OOM!
    // double min_cell_x = sph->box_size_x / 2000.0; // 確保 X 軸最多切 2000 格
    // double min_cell_y = sph->box_size_y / 2000.0;
    // double safe_min_h = fmax(min_cell_x, min_cell_y);
    // safe_min_h = fmax(safe_min_h, 1e-4); // 保留原本的絕對底線

    // max_h = max(max_h, safe_min_h);
    
    // 2. cell size
    sph->cell_size = max_h; 
    
    // 3. number of cells
    sph->num_cells_x = (int)ceil(sph->box_size_x / sph->cell_size);
    sph->num_cells_y = (int)ceil(sph->box_size_y / sph->cell_size);
    int new_total_cells = sph->num_cells_x * sph->num_cells_y;
    
    // 4. reallocate the length of the head array
    if (new_total_cells != sph->total_cells) {
        if (sph->head) free(sph->head);
        sph->head = (int *)malloc(new_total_cells * sizeof(int));
        sph->total_cells = new_total_cells;
    }
    
    // 5. initialize the head
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int c = 0; c < sph->total_cells; c++) {
        sph->head[c] = -1;
    }
    
    // 6. Create the Linked List
    // Do not parallelization!! Race Condition！
    // O(N), fast!
    for (int i = 0; i < sph->N; i++) {
        Particle *p = &sph->particles[i];
        
        // make sure it is inside
        double x = fmax(0.0, fmin(p->x, sph->box_size_x - 1e-9));
        double y = fmax(0.0, fmin(p->y, sph->box_size_y - 1e-9));
        
        // grid coordinate (cx, cy)
        int cx = (int)(x / sph->cell_size);
        int cy = (int)(y / sph->cell_size);
        
        // 1D array index
        int cell_index = cx + cy * sph->num_cells_x;
        
        // --- core Linked List insertion
        // 1. link the original head with i
        sph->next[i] = sph->head[cell_index]; 
        // 2. the head is now i
        sph->head[cell_index] = i; 
    }
}


void build_cell_list_3d(SPHSystem *sph) {
    double max_h = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(max:max_h)
#endif
    for (int i = 0; i < sph->N; i++) {
        if (sph->particles[i].h > max_h) {
            max_h = sph->particles[i].h;
        }
    }

    if (max_h < 1.0e-4) {
        max_h = 1.0e-4;
    }

    sph->cell_size = max_h;

    sph->num_cells_x = (int)ceil(sph->box_size_x / sph->cell_size);
    sph->num_cells_y = (int)ceil(sph->box_size_y / sph->cell_size);
    sph->num_cells_z = (int)ceil(sph->box_size_z / sph->cell_size);

    int new_total_cells =
        sph->num_cells_x * sph->num_cells_y * sph->num_cells_z;

    if (new_total_cells != sph->total_cells) {
        if (sph->head) {
            free(sph->head);
        }

        sph->head = (int *)malloc(new_total_cells * sizeof(int));

        if (sph->head == NULL) {
            fprintf(stderr, "Error: failed to allocate head in build_cell_list_3d.\n");
            exit(EXIT_FAILURE);
        }

        sph->total_cells = new_total_cells;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int c = 0; c < sph->total_cells; c++) {
        sph->head[c] = -1;
    }

    // Linked-list insertion: do not parallelize, otherwise race condition
    for (int i = 0; i < sph->N; i++) {
        Particle *p = &sph->particles[i];

        double x = fmax(0.0, fmin(p->x, sph->box_size_x - 1.0e-9));
        double y = fmax(0.0, fmin(p->y, sph->box_size_y - 1.0e-9));
        double z = fmax(0.0, fmin(p->z, sph->box_size_z - 1.0e-9));

        int cx = (int)(x / sph->cell_size);
        int cy = (int)(y / sph->cell_size);
        int cz = (int)(z / sph->cell_size);

        int cell_index = cx + cy * sph->num_cells_x + cz * sph->num_cells_x * sph->num_cells_y;

        sph->next[i] = sph->head[cell_index];
        sph->head[cell_index] = i;
    }
}
