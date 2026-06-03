#include "sph_all.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <chrono>
#endif

int main(int argc, char *argv[]) {
    printf("====================================\n");
    printf("   SPH Simulation Starting...             \n");
    printf("====================================\n");

    char output_folder[128];
    int custom_folder_set = 0;
    double mass = 0.001; 
    const char *mass_c = "0.001";
    double x = 15.0;
    const char *x_c = "15.0";
    double t_end = 5.0;
    const char *t_c = "5.0";
#ifndef __CUDACC__
    int num_threads = 0;
#endif
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
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            t_c = argv[i+1];
            t_end = atof(argv[i + 1]); 
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
#ifndef __CUDACC__
        // 新增：處理 -th 參數
        else if (strcmp(argv[i], "-th") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[i+1]);
            i++;
        }
#endif
    }

    if (custom_folder_set == 0){
        snprintf(output_folder, sizeof(output_folder), "sod_m%s_x%s_t%s", mass_c, x_c, t_c);
    }
    printf("Output will be saved to: %s\n", output_folder);

#ifdef __CUDACC__
    printf("Execution Mode: CUDA GPU Acceleration Enabled\n"); 
    int deviceId;
    cudaGetDevice(&deviceId);
    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, deviceId);
    printf("GPU Device: %s\n", props.name);
#elif defined(_OPENMP)
    if (num_threads > 0) {
        // 如果使用者有指定，就強制設定 OpenMP 使用該數量
        omp_set_num_threads(num_threads);
    }

    // 注意：在平行區塊「外」要使用 omp_get_max_threads()
    // 才能正確抓到「等一下平行區塊會用多少個 thread」
    int actual_threads = omp_get_max_threads();
    printf("Execution Mode: OpenMP Acceleration Enabled\n");
    printf("Number of Threads: %d\n", actual_threads);
#else
    if (num_threads > 0) {
        printf("Warning: '-th %d' ignored because OpenMP is not enabled during compilation.\n", num_threads);
    }
    printf("Execution Mode: Single Core (OpenMP disabled)\n");
#endif
    printf("====================================\n");

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

    // time_log.csv
    char log_filename[256];
    snprintf(log_filename, sizeof(log_filename), "%s/time_log.csv", output_folder);
    FILE *time_log = fopen(log_filename, "w");

    if (time_log != NULL) {
        // csv header
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
    // timing
#ifdef __CUDACC__
    auto start_time_gpu = std::chrono::high_resolution_clock::now();

#elif defined(_OPENMP)
    double start_time_omp=0.0, new_step_time=0.0, prev_step_time = 0.0;
    start_time_omp = omp_get_wtime();
    printf("Starting simulation with OpenMP acceleration...\n");
#else
    struct timeval run_t_start, run_t_init, run_t_end, run_t_new, run_t_old;
    gettimeofday(&run_t_start, NULL);
    printf("Starting simulation on Single Core...\n");
#endif

    // sod_2d_3 init
    init_sod_2d_3(&sph, x, y, mass);

    // Initialize physical parameters
    sph.gamma = 1.4;
    sph.cfl = 0.3;
    sph.alpha = 1.0;
    sph.beta = 2.0;
    sph.epsilon = 0.01;

    double eta = 2.3; 
    for (int i = 0; i < sph.N; i++) {
        Particle *p = &sph.particles[i];
        double target_rho = (p->x < 0.5 * x) ? 1.0 : 0.125; 

        // 直接套用 2D adaptive h 倒推公式
        p->h = eta * sqrt(p->mass / target_rho);
    }

    // Initial hydrodynamic computation before taking the first step
#ifdef __CUDACC__
    copy_particles_H2D(&sph); 

    compute_density_xreflective_yperiodic_celllist_gpu(&sph);
    compute_pressure_soundspeed_factor_gpu(&sph);
    compute_force_xreflective_yperiodic_celllist_gpu(&sph);
#else
    compute_density_xreflective_yperiodic_celllist(&sph);
    compute_pressure_soundspeed_factor(&sph);
    compute_force_xreflective_yperiodic_celllist(&sph);
#endif


    double t = 0.0;
    int step = 0;
    int output_step = 0;
    double dt_output = 0.01;
    double next_output_time = 0.0;


    double init_time = 0.0;
#ifdef __CUDACC__
    auto init_time_gpu = std::chrono::high_resolution_clock::now();
    auto prev_time_gpu = init_time_gpu;
#elif defined(_OPENMP)
    init_time = omp_get_wtime() - start_time_omp;
    prev_step_time = init_time;
#else
    gettimeofday(&run_t_init, NULL);
    init_time = (run_t_init.tv_sec - run_t_start.tv_sec) + (run_t_init.tv_usec - run_t_start.tv_usec) / 1000000.0;
    run_t_old = run_t_init;
#endif
    printf("Initialized in: %.3f seconds, start simulation\n", init_time);

    while (t < t_end + dt_output) {
        if (t >= next_output_time - 1e-9) {
            printf("output_time: %.4f\n", t);
            char filename[256];
            snprintf(filename, sizeof(filename), "%s/output_%04d.csv", output_folder, output_step);
#ifdef __CUDACC__
            copy_particles_D2H(&sph);
#endif
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
        if (t >= t_end - 1e-9)
            break;

        // integrate one step
#ifdef __CUDACC__
        double dt = step_leapfrog_kdk_xreflective_yperiodic_gpu(
                &sph, compute_timestep_signal_velocity_gpu, compute_force_xreflective_yperiodic_celllist_gpu);
#else
        double dt = step_leapfrog_kdk_xreflective_yperiodic(
                &sph, compute_timestep_signal_velocity, compute_force_xreflective_yperiodic_celllist);
#endif
        t += dt;
        step++;
        double step_time=0.0;
#ifdef __CUDACC__
        auto new_time_gpu = std::chrono::high_resolution_clock::now();
        step_time = std::chrono::duration<double>(new_time_gpu - prev_time_gpu).count();
        prev_time_gpu = new_time_gpu;
#elif defined(_OPENMP)
        new_step_time = omp_get_wtime();
        step_time = new_step_time - prev_step_time;
        prev_step_time = new_step_time;
#else
        gettimeofday(&run_t_new, NULL);
        step_time = (run_t_new.tv_sec - run_t_old.tv_sec) + (run_t_new.tv_usec - run_t_old.tv_usec) / 1000000.0;
        run_t_old = run_t_new;
#endif
        printf("sod t:%.10f, dt: %.10f, step simulation time: %f\n", t, dt, step_time);
    }

    double elapsed_time = 0.0;

#ifdef __CUDACC__
    elapsed_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time_gpu).count();
#elif defined(_OPENMP)
    elapsed_time = omp_get_wtime() - start_time_omp;
#else
    gettimeofday(&run_t_end, NULL);
    elapsed_time = (run_t_end.tv_sec - run_t_start.tv_sec) + (run_t_end.tv_usec - run_t_start.tv_usec) / 1000000.0;
#endif

    printf("====================================\n");
    printf("Simulation finished!\n");
#ifdef __CUDACC__
    printf("GPU Device: %s\n", props.name);
#endif
    printf("Total Execution Time: %.3f seconds\n", elapsed_time);
    printf("Initialize Time: %.3f seconds\n", init_time);
    printf("Simulation Time: %.3f seconds\n", elapsed_time - init_time);
    printf("Mass of particles: %f\n", mass);
    printf("Number of particles: %d\n", sph.N);
    printf("Output saved to: %s\n", output_folder);
    printf("Total Outputfile numbers: %d\n", output_step);
    printf("====================================\n");

    // 3. 迴圈結束後，安全關閉檔案
    if (time_log != NULL)fclose(time_log);
    free_sph_system(&sph);
    return 0;
}
