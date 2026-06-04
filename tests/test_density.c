#include "sph_all.h"

int main()
{
    printf("========== SPH Density Test ==========\n");

    SPHSystem sph;

    allocate_sph_system(&sph, 2);

    sph.particles[0].id = 0;
    sph.particles[0].x = 0.0;
    sph.particles[0].y = 0.0;
    sph.particles[0].mass = 1.0;
    sph.particles[0].h = 1.0;
    sph.particles[0].rho = 0.0;

    sph.particles[1].id = 1;
    sph.particles[1].x = 0.5;
    sph.particles[1].y = 0.0;
    sph.particles[1].mass = 1.0;
    sph.particles[1].h = 1.0;
    sph.particles[1].rho = 0.0;

    compute_density(&sph);

    printf("Particle 0 Density: %f\n", sph.particles[0].rho);
    printf("Particle 1 Density: %f\n", sph.particles[1].rho);

    printf("Note: densities should be the same.\n");

    free_sph_system(&sph);

    printf("========== Test Finished ==========\n");

    return 0;
}
