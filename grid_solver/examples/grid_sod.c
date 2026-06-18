#include "../include/grid_init.h"
#include "../include/grid_io.h"
#include "../include/grid_system.h"
#include "../include/hydro.h"
#include <sys/time.h>

#include <string.h>
#include <errno.h>
#include <sys/stat.h>

int main(int argc, char *argv[]) {
  printf("====================================\n");
  printf("   1D Grid Sod Shock Tube         \n");
  printf("====================================\n");

  int nx = 15000;
  char output_folder[128] = "grid_sod_output";

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      nx = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      strncpy(output_folder, argv[i + 1], sizeof(output_folder) - 1);
      output_folder[sizeof(output_folder) - 1] = '\0';
      i++;
    }
  }

  if (mkdir(output_folder, 0777) == -1) {
    if (errno != EEXIST) {
      printf("Error: Could not create directory %s\n", output_folder);
      return 1;
    }
  }

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
      char filename[256];
      sprintf(filename, "%s/output_%04d.csv", output_folder, output_step++);
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
  char filename[256];
  sprintf(filename, "%s/output_%04d.csv", output_folder, output_step);
  write_csv(&grid, filename);

  gettimeofday(&end, NULL);
  double elapsed =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
  printf("Simulation completed in %.3f seconds, %d steps.\n", elapsed, step);

  free_grid(&grid);
  return 0;
}
