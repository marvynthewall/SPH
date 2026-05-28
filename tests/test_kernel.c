#include "../include/sph_system.h"

int main() {
  printf("========== SPH Kernel Test ==========\n");

  double h = 1.0;
  double W, dWdr, dWdh;
  printf("\nTesting Kernel (h = %.1f)\n", h);

  double sigma = 40.0 / (7.0 * M_PI * h * h);
  // 測試 r = 0 (自己對自己)
  cubic_spline_kernel_2d(0.0, h, &W, &dWdr, &dWdh);
  printf("At r = 0.0: W = %f (Expected max value, sigma = %f), dWdr = %f\n", W, sigma, dWdr);

  // 測試 r = 0.5 (在影響範圍內)
  cubic_spline_kernel_2d(0.5, h, &W, &dWdr, &dWdh);
  printf("At r = 0.5: W = %f, dWdr = %f\n", W, dWdr);

  // 測試 r = 1.5 (超過影響範圍 h，權重應該要是 0)
  cubic_spline_kernel_2d(1.5, h, &W, &dWdr, &dWdh);
  printf("At r = 1.5: W = %f (Expected 0.0), dWdr = %f\n", W, dWdr);

  printf("\n================ Test Finished ================\n");
  return 0;
}
