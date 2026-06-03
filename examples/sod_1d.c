#include "density.h"
#include "force.h"
#include "init.h"
#include "integrator.h"
#include "io.h"
#include "sph_system.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <errno.h>
#include <sys/stat.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
  printf("====================================\n");
  printf("   1D SPH Simulation Starting...    \n");
  printf("====================================\n");

  char output_folder[128];
  int custom_folder_set = 0;
  double mass = 0.001;
  const char *mass_c = "0.001";
  double x = 15.0;
  const char *x_c = "15.0";
  double t_end = 5.0;
  const char *t_c = "5.0";
  int num_threads = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
      mass_c = argv[i + 1];
      mass = atof(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-x") == 0 && i + 1 < argc) {
      x_c = argv[i + 1];
      x = atof(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
      t_c = argv[i + 1];
      t_end = atof(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      strncpy(output_folder, argv[i + 1], sizeof(output_folder) - 1);
      output_folder[sizeof(output_folder) - 1] = '\0';
      custom_folder_set = 1;
      i++;
    } else if (strcmp(argv[i], "-th") == 0 && i + 1 < argc) {
      num_threads = atoi(argv[i + 1]);
      i++;
    }
  }

  if (custom_folder_set == 0) {
    snprintf(output_folder, sizeof(output_folder), "sod_1d_m%s_x%s_t%s", mass_c,
             x_c, t_c);
  }
  printf("Output will be saved to: %s\n", output_folder);

#ifdef _OPENMP
  if (num_threads > 0) {
    omp_set_num_threads(num_threads);
  }
  int actual_threads = omp_get_max_threads();
  printf("Execution Mode: OpenMP Acceleration Enabled\n");
  printf("Number of Threads: %d\n", actual_threads);
#else
  printf("Execution Mode: Single Core (OpenMP disabled)\n");
#endif
  printf("====================================\n");

  if (mkdir(output_folder, 0777) == -1) {
    if (errno != EEXIST) {
      printf("Error: Could not create directory %s\n", output_folder);
      return 1;
    }
  }

  SPHSystem sph;

  char log_filename[256];
  snprintf(log_filename, sizeof(log_filename), "%s/time_log.csv",
           output_folder);
  FILE *time_log = fopen(log_filename, "w");

  if (time_log != NULL) {
    fprintf(time_log, "frame,t\n");
  } else {
    printf("Warning: Could not open %s for writing.\n", log_filename);
  }

#ifdef _OPENMP
  double start_time_omp = 0.0, new_step_time = 0.0, prev_step_time = 0.0;
  start_time_omp = omp_get_wtime();
  printf("Starting simulation with OpenMP acceleration...\n");
#else
  struct timeval run_t_start, run_t_init, run_t_end, run_t_new, run_t_old;
  gettimeofday(&run_t_start, NULL);
  printf("Starting simulation on Single Core...\n");
#endif

  // sod_1d init
  init_sod_1d(&sph, x, mass);

  sph.gamma = 1.4;
  sph.cfl = 0.3;
  sph.alpha = 1.0;
  sph.beta = 2.0;
  sph.epsilon = 0.01;
  sph.dim = 1; // Used in compute_pressure_soundspeed_factor

  // Initial hydrodynamic computation before taking the first step
  compute_density_1d_xreflective(&sph);
  compute_pressure_soundspeed_factor(&sph);
  compute_force_1d_xreflective(&sph);

  double t = 0.0;
  int step = 0;
  int output_step = 0;
  double dt_output = 0.01;
  double next_output_time = 0.0;

  double init_time = 0.0;
#ifdef _OPENMP
  init_time = omp_get_wtime() - start_time_omp;
#else
  gettimeofday(&run_t_init, NULL);
  init_time = (run_t_init.tv_sec - run_t_start.tv_sec) +
              (run_t_init.tv_usec - run_t_start.tv_usec) / 1000000.0;
  run_t_old = run_t_init;
#endif
  printf("Initialized in: %.3f seconds, start simulation\n", init_time);

  while (t < t_end + dt_output) {
    if (t >= next_output_time - 1e-9) {
      char filename[256];
      snprintf(filename, sizeof(filename), "%s/output_%04d.csv", output_folder,
               output_step);
      write_csv(&sph, filename);
      printf("Step %d | Time: %.4f | dt: %.6f | Output: %s\n", step, t, sph.dt,
             filename);
      if (time_log != NULL) {
        fprintf(time_log, "%d,%.10f\n", output_step, t);
        fflush(time_log);
      }
      output_step++;
      next_output_time += dt_output;
    }
    if (t >= t_end - 1e-9)
      break;

    double dt = step_leapfrog_kdk_1d_xreflective(
        &sph, compute_timestep_signal_velocity, compute_force_1d_xreflective);
    t += dt;
    step++;
    double step_time;
#ifdef _OPENMP
    new_step_time = omp_get_wtime();
    step_time = new_step_time - prev_step_time;
    prev_step_time = new_step_time;
#else
    gettimeofday(&run_t_new, NULL);
    step_time = (run_t_new.tv_sec - run_t_old.tv_sec) +
                (run_t_new.tv_usec - run_t_old.tv_usec) / 1000000.0;
    run_t_old = run_t_new;
#endif
  }

  double elapsed_time = 0.0;

#ifdef _OPENMP
  elapsed_time = omp_get_wtime() - start_time_omp;
#else
  gettimeofday(&run_t_end, NULL);
  elapsed_time = (run_t_end.tv_sec - run_t_start.tv_sec) +
                 (run_t_end.tv_usec - run_t_start.tv_usec) / 1000000.0;
#endif

  printf("====================================\n");
  printf("Simulation finished!\n");
  printf("Total Execution Time: %.3f seconds\n", elapsed_time);
  printf("Initialize Time: %.3f seconds\n", init_time);
  printf("Simulation Time: %.3f seconds\n", elapsed_time - init_time);
  printf("Total Outputfile numbers: %d\n", output_step);
  printf("Mass of particles: %f\n", mass);
  printf("====================================\n");

  if (time_log != NULL)
    fclose(time_log);
  free_sph_system(&sph);
  return 0;
}
