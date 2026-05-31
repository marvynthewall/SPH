#include "sph_system.h"

int main(int argc, char *argv[]) {
  printf("====================================\n");
  printf("   KH Instability 模擬程式啟動中...             \n");
  printf("====================================\n");

  int nx = 128;
  int ny = 128;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      nx = atoi(argv[i + 1]);
      ny = atoi(argv[i + 1]);
      i++;
    }
  }

  SPHSystem2D sph;
  
  // Initialize Kelvin-Helmholtz instability setup
  init_KH(&sph, nx, ny);

  // Initialize physical parameters
  sph.cfl = 0.3;
  sph.alpha = 1.0;
  sph.beta = 2.0;
  sph.epsilon = 0.01;

  // Initial hydrodynamic computation before taking the first step
  compute_density_xperiodic_yperiodic(&sph);
  compute_pressure_soundspeed_factor(&sph);
  compute_force_xperiodic_yperiodic(&sph);

  double t = 0.0;
  double t_end = 2.0;
  int step = 0;
  int output_step = 0;
  double dt_output = 0.01;
  double next_output_time = 0.0;

  printf("\n初始化完成,開始模擬...\n");

  while (t < t_end) {
    if (t >= next_output_time) {
      printf("output_time: %.4f\n", t);
      char filename[256];
      sprintf(filename, "output_%04d.csv", output_step);
      write_csv(&sph, filename);
      printf("Step %d | Time: %.4f | dt: %.6f | Output: %s\n", step, t, sph.dt,
             filename);
      output_step++;
      next_output_time += dt_output;
    }

    // integrate one step using fully periodic boundary conditions
    double dt = step_leapfrog_kdk_xperiodic_yperiodic(
        &sph, compute_timestep_signal_velocity, compute_force_xperiodic_yperiodic);
    t += dt;
    step++;
  }

  printf("\n模擬完成 ...\n");
  printf("總共輸出檔案數量: %d\n", output_step);
  printf("設定網格解析度: %dx%d\n", nx, ny);

  free_sph_system(&sph);
  return 0;
}
