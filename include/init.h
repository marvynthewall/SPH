#ifndef INIT_H
#define INIT_H

#include "sph_system.h"

/**
 * @brief Initializes particles in a 2D square region.
 * @param particles Array of particles to initialize.
 * @param num_particles Number of particles to create.
 * @note Particles are initialized in a grid pattern with random perturbations.
 */
void initialize_particles(Particle *particles, int num_particles);

#endif