#include "../include/kernel.h"

void cubic_spline_kernel(double r, double h, double *W, double *dWdr, double *dWdh)
{
  // q is the dimensionless distance (r/h)
  double q = r / h;

  // 2D normalization constant sigma = 40 / (7 * pi)
  // Add 1/h^2 to ensure the units and integral are correct
  double norm = 40.0 / (7.0 * M_PI * h * h);

  // Pre-calculate the constant term for the derivative, because
  // dW/dr = (dW/dq) * (dq/dr) = (dW/dq) * (1/h)
  double norm_dW = norm / h;

  // (d_/d_) means partial derivative
  // let k = norm
  // W = k(h) * P(q)
  // dW/dr = (dW/dq) * (1/h) = (k/h) * dP/dq
  // k = norm = a / h^2 => dk/dh = (-2/h) * k
  // q = r / h          => dq/dh = -r / h^2
  // dW/dh  = (dk/dh) * P(q) + k * (dP/dq) * (dq/dr)
  //        = (-2/h) * W + k * (dP/dq) * (-r / h^2)
  //        = (-2/h) * W - (r/h) * (k/h) * (dP/dq)
  //        = (-2/h) * W - q * (dW/dr)

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
    // Beyond the influence range (q > 1 or r > h), the weight is 0
    *W = 0.0;
    *dWdr = 0.0;
    *dWdh = 0.0;
  }
}

void cubic_spline_kernel_3d(double r, double h, double *W, double *dWdr, double *dWdh)
{
  // q is the dimensionless distance (r/h)
  double q = r / h;

  // 3D normalization constant sigma = 8 / pi (Eq.2)
  // Multiply by 1/h^3 to ensure correct 3D volumetric units and integral
  double h3 = h * h * h;
  double norm = 8.0 / (M_PI * h3);

  // Pre-calculate the constant term for the radial derivative
  // dW/dr = (dW/dq) * (dq/dr) = (dW/dq) * (1/h)
  double norm_dW = norm / h;

  // Piecewise polynomial calculation
  if (q >= 0.0 && q <= 0.5) {
    *W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    *dWdr = norm_dW * (-12.0 * q + 18.0 * q * q);

    // 3D derivative w.r.t h uses -3/h (volumetric scaling) instead of -2/h
    // (areal)
    *dWdh = (-3.0 / h) * (*W) - q * (*dWdr);

  } else if (q > 0.5 && q <= 1.0) {
    double diff = 1.0 - q;
    *W = norm * 2.0 * diff * diff * diff;
    *dWdr = norm_dW * (-6.0 * diff * diff);

    *dWdh = (-3.0 / h) * (*W) - q * (*dWdr);

  } else {
    // Beyond the influence range (q > 1 or r > h), the weight is 0
    *W = 0.0;
    *dWdr = 0.0;
    *dWdh = 0.0;
  }
}

void cubic_spline_kernel_1d(double r, double h, double *W, double *dWdr, double *dWdh)
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
