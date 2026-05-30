#include "sph_system.h"

int main(int argc, char *argv[]) {
  printf("====================================\n");
  printf("   SPH 模擬程式啟動中...             \n");
  printf("====================================\n");
  // 1. 設定預設值
  const char *output_filename = "output_0000.csv";
  double mass = 0.005; // 預設質量

  // 2. 解析命令列參數
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      output_filename = argv[i + 1];
      i++; // 跳過已讀取的檔名
    } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
      mass = atof(argv[i + 1]); // 將字串轉換為 double
      i++;                      // 跳過已讀取的數值
    }
  }

  SPHSystem2D sph;
  double x = 1.0;
  double y = 1.0;

  // 如果你的 init_sod_2d_2 需要用到 mass，記得要把 mass 傳進去，
  // 例如：init_sod_2d_2(&sph, x, y, mass);
  // 這裡先維持你原本的呼叫方式：
  // init_sod_2d_2(&sph, x, y, mass);

  // uniform init
  /*
  allocate_sph_system(&sph, 256);
  sph.box_size_x = 1.0;
  sph.box_size_y = 1.0;
  init_uniform_box(&sph, 16, 16);
  */

  // sod_2d_3 init
  init_sod_2d_3(&sph, x, y, mass);

  // Initialize physical parameters
  sph.cfl = 0.3;
  sph.alpha = 1.0;
  sph.beta = 2.0;
  sph.epsilon = 0.01;

  // Initial hydrodynamic computation before taking the first step
  compute_density_xreflective_yperiodic(&sph);
  compute_pressure_soundspeed_factor(&sph);
  compute_force_xreflective_yperiodic(&sph);

  double t = 0.0;
  double t_end = 0.5;
  int step = 0;
  int output_step = 0;
  double dt_output = 0.01;
  double next_output_time = 0.0;

  printf("\n初始化完成，開始模擬...\n");

  while (t < t_end) {
    if (t >= next_output_time) {
      char filename[256];
      sprintf(filename, "output_%04d.csv", output_step);
      write_csv(&sph, filename);
      printf("Step %d | Time: %.4f | dt: %.6f | Output: %s\n", step, t, sph.dt,
             filename);
      output_step++;
      next_output_time += dt_output;
    }

    // integrate one step
    double dt = step_leapfrog_kdk_xreflective_yperiodic(
        &sph, compute_timestep, compute_force_xreflective_yperiodic);
    t += dt;
    step++;
  }

  printf("\n模擬完成 ...\n");
  printf("總共輸出檔案數量: %d\n", output_step);
  printf("設定粒子質量: %f\n", mass);

  free_sph_system(&sph);
  return 0;
}
