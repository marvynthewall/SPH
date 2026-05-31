#include "sph_system.h"
#include "density.h"
#include "init.h"


int main(int argc, char *argv[]) {
    SPHSystem2D sph;

    int nx = 20; 
    int ny = 20;
    int N  = nx * ny;
    double eta = 2.3;

    allocate_sph_system(&sph, N);
    init_uniform_box(&sph, nx, ny);

    update_adaptive_h_2d(&sph, 1000, 1.0e-4, eta, compute_density_xreflective_yperiodic);
    // check_adaptive_h(&sph, 2.3, 1.0e-4);
    free_sph_system(&sph);
    
    return 0;
}