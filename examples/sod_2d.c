#include "sph_system.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

int main(int argc, char *argv[]) {
  printf("====================================\n");
  printf("   SPH 模擬程式啟動中...             \n");
  printf("====================================\n");

  char output_folder[128];
  int custom_folder_set = 0;
  double mass = 0.001; 
  char *mass_c = "0.001";
  double x = 5.0;
  char *x_c = "5.0";

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
      mass_c = argv[i+1];
      mass = atof(argv[i + 1]); 
      i++;                      
    }
    else if (strcmp(argv[i], "-x") == 0 && i + 1 < argc) {
      x_c = argv[i+1];
      x = atof(argv[i + 1]); 
      i++;                      
    }
    // 新增 -o 參數的判斷
    else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      // 使用 strncpy 安全地複製字串，預留 1 byte 給結尾符號 '\0'
      strncpy(output_folder, argv[i+1], sizeof(output_folder) - 1);
      output_folder[sizeof(output_folder) - 1] = '\0'; // 強制加上字串結尾符號

      custom_folder_set = 1; // 標記使用者有指定資料夾
      i++;
    }
  }
  
  if (custom_folder_set == 0){
    snprintf(output_folder, sizeof(output_folder), "sod_default_folder");
  }
  snprintf(output_folder, sizeof(output_folder), "sod_m%s_x%s", mass_c, x_c);
  printf("Output will be saved to: %s\n", output_folder);

  // 嘗試建立資料夾，0777 是基本的讀寫權限設定
  if (mkdir(output_folder, 0777) == -1) {
      // 如果建立失敗，且原因不是「資料夾已經存在 (EEXIST)」
      if (errno != EEXIST) {
          printf("Error: Could not create directory %s\n", output_folder);
          return 1; // 終止程式
      }
  }

  SPHSystem sph;
  double y = 1.0;

  // 1. 在進入迴圈前開啟 time_log.csv
  char log_filename[256];
  snprintf(log_filename, sizeof(log_filename), "%s/time_log.csv", output_folder);
  FILE *time_log = fopen(log_filename, "w");

  if (time_log != NULL) {
      // 寫入 CSV 標頭
      fprintf(time_log, "frame,t\n");
  } else {
      printf("Warning: Could not open %s for writing.\n", log_filename);
  }


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
  sph.gamma = 1.4;
  sph.cfl = 0.3;
  sph.alpha = 1.0;
  sph.beta = 2.0;
  sph.epsilon = 0.01;

  // Initial hydrodynamic computation before taking the first step
  compute_density_xreflective_yperiodic(&sph);
  compute_pressure_soundspeed_factor(&sph);
  compute_force_xreflective_yperiodic(&sph);

  double t = 0.0;
  double t_end = 5.0;
  int step = 0;
  int output_step = 0;
  double dt_output = 0.01;
  double next_output_time = 0.0;

  printf("\n初始化完成,開始模擬...\n");

  while (t < t_end) {
    if (t >= next_output_time) {
      printf("output_time: %.4f\n", t);
      char filename[256];
      snprintf(filename, sizeof(filename), "%s/output_%04d.csv", output_folder, output_step);
      write_csv(&sph, filename);
      printf("Step %d | Time: %.4f | dt: %.6f | Output: %s\n", step, t, sph.dt,
             filename);
      // 2. 將當下的 frame 與精確時間 t 寫入 log
      if (time_log != NULL) {
          fprintf(time_log, "%d,%.10f\n", output_step, t);
          // 加上 fflush 強制寫入硬碟，確保模擬如果中途中斷，已經跑完的時間紀錄不會遺失
          fflush(time_log);
      }
      output_step++;
      next_output_time += dt_output;

    }

    // integrate one step
    // double dt = step_leapfrog_kdk_xreflective_yperiodic(
    //     &sph, compute_timestep, compute_force_xreflective_yperiodic);
    double dt = step_leapfrog_kdk_xreflective_yperiodic(
        &sph, compute_timestep_signal_velocity, compute_force_xreflective_yperiodic);
    t += dt;
    printf("sod t:%.10f, dt: %.10f\n", t, dt);
    step++;
  }

  printf("\n模擬完成 ...\n");
  printf("總共輸出檔案數量: %d\n", output_step);
  printf("設定粒子質量: %f\n", mass);

  // 3. 迴圈結束後，安全關閉檔案
  if (time_log != NULL)fclose(time_log);
  free_sph_system(&sph);
  return 0;
}
