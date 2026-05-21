#ifndef KERNEL_H
#define KERNEL_H

#include <math.h>

/*
 * Cubic Spline Kernel for 2D
 *
 * Physical Meaning:
 * This is the classic weighting function in SPH.
 * It discretizes a continuous fluid field to individual particles.
 * Neighbors closer to the center particle contribute more to its physical
 * quantities (like mass). When the distance exceeds h, the influence drops to 0
 * (Compact support). Note: While some literature defines the cutoff at 2h, in
 * our formulation where q = r/h and the kernel vanishes for q > 1, the cutoff
 * distance is exactly h.
 *
 * Computational Technique:
 * To avoid repeated calculations, we define the dimensionless distance as q = r
 * / h. We return both W (kernel value) and dWdr (derivative with respect to r)
 * simultaneously, as the derivative is needed for force calculations (e.g.,
 * pressure gradient).
 *
 * Parameters:
 * r: Actual distance between two particles |r_i - r_j|
 * h: Smoothing length, determining the influence range of the particle
 * W: Calculated kernel value (output)
 * dWdr: Calculated derivative of the kernel with respect to distance r (output)
 */
void cubic_spline_kernel_2d(double r, double h, double *W, double *dWdr);

#endif