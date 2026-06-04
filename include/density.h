#pragma once
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
void compute_density(SPHSystem *sph);
void compute_density_3d(SPHSystem *sph);
void compute_density_1d_xreflective(SPHSystem *sph);
void compute_density_xreflective_yperiodic(SPHSystem *sph);
void compute_density_xperiodic_yperiodic(SPHSystem *sph);
void compute_density_xreflective_yperiodic_celllist(SPHSystem *sph);
void compute_density_xreflective_yzperiodic_3d(SPHSystem *sph);
void compute_density_xreflective_yzperiodic_celllist_3d(SPHSystem *sph);

void update_adaptive_h(SPHSystem *sph, int max_iter, double tol, double eta,
                       void (*compute_density_fn)(SPHSystem *));
void update_adaptive_h_3d(SPHSystem *sph, int max_iter, double tol, double eta,
                          void (*compute_density_fn)(SPHSystem *));
void check_adaptive_h(SPHSystem *sph, double eta, double tol);

// GPU

void update_adaptive_h_gpu(SPHSystem *sph, int max_iter, double tol, double eta,
                       void (*compute_density_fn)(SPHSystem *));
void compute_density_xreflective_yperiodic_celllist_gpu(SPHSystem *sph);
