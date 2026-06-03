#pragma once
#include "sph_system.h"

#ifdef __CUDACC__
    #define SPH_MATH_INLINE __host__ __device__ __forceinline__
#else
    #define SPH_MATH_INLINE static inline
#endif
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
// void cubic_spline_kernel_1d(double r, double h, double *W, double *dWdr, double *dWdh);
// void cubic_spline_kernel(double r, double h, double *W, double *dWdr, double *dWdh);
// void cubic_spline_kernel_3d(double r, double h, double *W, double *dWdr, double *dWdh);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SPH_MATH_INLINE void cubic_spline_kernel_1d(double r, double h, double *W, double *dWdr, double *dWdh)
{
  double q = r / h;
  double norm = 4.0 / (3.0 * h);
  double norm_dW = norm / h;

  if (q >= 0.0 && q <= 0.5) {
    *W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    *dWdr = norm_dW * (-12.0 * q + 18.0 * q * q);
    *dWdh = (-1.0 / h) * (*W) - q * (*dWdr);
  } else if (q > 0.5 && q <= 1.0) {
    double diff = 1.0 - q;
    *W = norm * 2.0 * diff * diff * diff;
    *dWdr = norm_dW * (-6.0 * diff * diff);
    *dWdh = (-1.0 / h) * (*W) - q * (*dWdr);
  } else {
    *W = 0.0;
    *dWdr = 0.0;
    *dWdh = 0.0;
  }
}

SPH_MATH_INLINE void cubic_spline_kernel(double r, double h, double *W, double *dWdr, double *dWdh)
{
  // q is the dimensionless distance (r/h)
  double q = r / h;

  // 2D normalization constant sigma = 40 / (7 * pi)
  // Add 1/h^2 to ensure the units and integral are correct
  double norm = 40.0 / (7.0 * M_PI * h * h);

  double norm_dW = norm / h;

  // Piecewise polynomial calculation
  if (q >= 0.0 && q <= 0.5) {
    *W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    *dWdr = norm_dW * (-12.0 * q + 18.0 * q * q);
    *dWdh = (-2.0 / h) * (*W) - q * (*dWdr);
  } else if (q > 0.5 && q <= 1.0) {
    double diff = 1.0 - q;
    *W = norm * 2.0 * diff * diff * diff;
    *dWdr = norm_dW * (-6.0 * diff * diff);
    *dWdh = (-2.0 / h) * (*W) - q * (*dWdr);
  } else {
    *W = 0.0;
    *dWdr = 0.0;
    *dWdh = 0.0;
  }
}

SPH_MATH_INLINE void cubic_spline_kernel_3d(double r, double h, double *W, double *dWdr, double *dWdh)
{
  // q is the dimensionless distance (r/h)
  double q = r / h;

  // 3D normalization constant sigma = 8 / pi (Eq.2)
  double h3 = h * h * h;
  double norm = 8.0 / (M_PI * h3);
  double norm_dW = norm / h;

  // Piecewise polynomial calculation
  if (q >= 0.0 && q <= 0.5) {
    *W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    *dWdr = norm_dW * (-12.0 * q + 18.0 * q * q);
    *dWdh = (-3.0 / h) * (*W) - q * (*dWdr);
  } else if (q > 0.5 && q <= 1.0) {
    double diff = 1.0 - q;
    *W = norm * 2.0 * diff * diff * diff;
    *dWdr = norm_dW * (-6.0 * diff * diff);
    *dWdh = (-3.0 / h) * (*W) - q * (*dWdr);
  } else {
    *W = 0.0;
    *dWdr = 0.0;
    *dWdh = 0.0;
  }
}
