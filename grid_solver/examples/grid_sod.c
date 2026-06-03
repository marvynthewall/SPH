#include "../include/grid_init.h"
#include "../include/grid_io.h"
#include "../include/grid_system.h"
#include "../include/hydro.h"
#include <sys/time.h>

int main(void) {
  printf("====================================\n");
  printf("   1D Grid Sod Shock Tube         \n");
  printf("====================================\n");

  int nx = 474; // Resolution matching 2D SPH dx = 0.0316
  int ny = 1;
  int ng = 1; // 1 ghost cell for 1st order Godunov

  double xmin = 0.0, xmax = 15.0;
  double ymin = 0.0, ymax = 1.0; // dummy for 1D
  double gamma = 1.4;
  double cfl = 0.3;
  double t_end = 5.0;

  GridSystem grid;
  allocate_grid(&grid, nx, ny, ng, xmin, xmax, ymin, ymax, gamma, cfl, t_end);

  init_sod_1d(&grid);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  int step = 0;
  int output_step = 0;
  double output_dt = 0.01; // Output every 0.01 time units
  double next_output_t = 0.0;

  while (grid.t < grid.t_end) {
    if (grid.t >= next_output_t) {
      char filename[128];
      sprintf(filename, "grid_sod_2d_output/output_%04d.csv", output_step++);
      write_csv(&grid, filename);
      next_output_t += output_dt;
    }

    apply_boundary_conditions(&grid, 0); // 0 = Sod
    calculate_dt(&grid);
    hydro_step(&grid);

    grid.t += grid.dt;
    step++;

    if (step % 100 == 0) {
      printf("Step %d, t = %.4f, dt = %.4e\n", step, grid.t, grid.dt);
    }
  }

  // Final output
  char filename[128];
  sprintf(filename, "grid_sod_2d_output/output_%04d.csv", output_step);
  write_csv(&grid, filename);

  gettimeofday(&end, NULL);
  double elapsed =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
  printf("Simulation completed in %.3f seconds, %d steps.\n", elapsed, step);

  free_grid(&grid);
  return 0;
}
