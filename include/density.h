#ifndef DENSITY_H
#define DENSITY_H


#include "sph_system.h"
/*
 * Calculate the density of all particles
 *
 * Physical Meaning:
 * According to SPH theory, density is defined as: rho_i = sum_j ( m_j * W_ij )
 * This means the density of particle i is calculated by summing the masses of
 * neighboring particles (distance less than h) multiplied by the kernel
 * function weights. Note that a particle also contributes to its own density
 * (i.e., when j = i, using W(0, h)).
 *
 * Computational Technique:
 * If Neighbor List or Cell-Linked List is used for acceleration,
 * then each particle only needs to search within its own cell or adjacent
 * cells, significantly reducing the time complexity from O(N^2).
 *
 * Parameters:
 * particles: Array of all particles in the system
 * num_particles: Total number of particles
 */
void compute_density(SPHSystem2D *sph);
void compute_density_xreflective_yperiodic(SPHSystem2D *sph);

void update_adaptive_h_2d(SPHSystem2D *sph, int max_iter, double tol, double eta, void (*compute_density_fn)(SPHSystem2D *));
void update_adaptive_h_3d(SPHSystem2D *sph, int max_iter, double tol, double eta, void (*compute_density_fn)(SPHSystem2D *));
void check_adaptive_h(SPHSystem2D *sph, double eta, double tol);
#endif