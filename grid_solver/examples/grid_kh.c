#include "../include/grid_init.h"
#include "../include/grid_io.h"
#include "../include/grid_system.h"
#include "../include/hydro.h"
#include <sys/time.h>

int main(void) {
  printf("====================================\n");
  printf("   2D Grid Kelvin-Helmholtz       \n");
  printf("====================================\n");

  int nx = 128; // Resolution
  int ny = 128;
  int ng = 1; // 1 ghost cell for 1st order Godunov

  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double gamma = 5.0 / 3.0; // Monatomic gas for KH
  double cfl = 0.3;
  double t_end = 2.0; // typical KH simulation time

  GridSystem grid;
  allocate_grid(&grid, nx, ny, ng, xmin, xmax, ymin, ymax, gamma, cfl, t_end);

  init_kh_2d(&grid);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  int step = 0;
  int output_step = 0;
  double output_dt = 0.01;
  double next_output_t = 0.0;

  while (grid.t < grid.t_end) {
    if (grid.t >= next_output_t) {
      char filename[128];
      sprintf(filename, "grid_kh_2d_output/output_%04d.csv", output_step++);
      write_csv(&grid, filename);
      next_output_t += output_dt;
      printf("Output %d written at t = %.4f\n", output_step - 1, grid.t);
    }

    apply_boundary_conditions(&grid, 1); // 1 = KH
    calculate_dt(&grid);
    hydro_step(&grid);

    grid.t += grid.dt;
    step++;
  }

  // Final output
  char filename[128];
  sprintf(filename, "grid_kh_2d_output/output_%04d.csv", output_step);
  write_csv(&grid, filename);

  gettimeofday(&end, NULL);
  double elapsed =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
  printf("Simulation completed in %.3f seconds, %d steps.\n", elapsed, step);

  free_grid(&grid);
  return 0;
}
