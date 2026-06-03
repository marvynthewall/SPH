#include "sph_system.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
    printf("====================================\n");
    printf("   3D SPH Sod 模擬程式啟動中...      \n");
    printf("====================================\n");

    char output_folder[128];
    int custom_folder_set = 0;

    double mass = 0.001;
    char *mass_c = "0.001";

    double x = 15.0;
    const char *x_c = "15.0";

    double y = 1.0;
    char *y_c = "1.0";

    double z = 1.0;
    char *z_c = "1.0";

    double t_end = 20.0;
    char *t_c = "5.0";

    int num_threads = 0;

    for (int i = 1; i < argc; i++) {

        if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            mass_c = argv[i + 1];
            mass = atof(argv[i + 1]);
            i++;
        }

        else if (strcmp(argv[i], "-x") == 0 && i + 1 < argc) {
            x_c = argv[i + 1];
            x = atof(argv[i + 1]);
            i++;
        }

        else if (strcmp(argv[i], "-y") == 0 && i + 1 < argc) {
            y_c = argv[i + 1];
            y = atof(argv[i + 1]);
            i++;
        }

        else if (strcmp(argv[i], "-z") == 0 && i + 1 < argc) {
            z_c = argv[i + 1];
            z = atof(argv[i + 1]);
            i++;
        }

        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            t_c = argv[i + 1];
            t_end = atof(argv[i + 1]);
            i++;
        }

        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            strncpy(output_folder, argv[i + 1], sizeof(output_folder) - 1);
            output_folder[sizeof(output_folder) - 1] = '\0';
            custom_folder_set = 1;
            i++;
        }

        else if (strcmp(argv[i], "-th") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[i + 1]);
            i++;
        }
    }

    if (custom_folder_set == 0) {
        snprintf(output_folder,
                 sizeof(output_folder),
                 "sod3d_m%s_x%s_y%s_z%s_t%s",
                 mass_c,
                 x_c,
                 y_c,
                 z_c,
                 t_c);
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
    if (num_threads > 0) {
        printf("Warning: '-th %d' ignored because OpenMP is not enabled.\n",
               num_threads);
    }

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
    snprintf(log_filename, sizeof(log_filename), "%s/time_log.csv", output_folder);

    FILE *time_log = fopen(log_filename, "w");

    if (time_log != NULL) {
        fprintf(time_log, "frame,t\n");
    } else {
        printf("Warning: Could not open %s for writing.\n", log_filename);
    }

    init_sod_3d_3(&sph, x, y, z, mass);

    sph.gamma = 1.4;
    sph.cfl = 0.1;
    sph.alpha = 0.5;
    sph.beta = 0.3;
    sph.epsilon = 0.01;

    compute_density_xreflective_yzperiodic_3d(&sph);
    compute_pressure_soundspeed_factor(&sph);
    compute_force_xreflective_yperiodic_zperiodic_3d(&sph);

    double t = 0.0;
    int step = 0;
    int output_step = 0;

    double dt_output = 1.0;
    double next_output_time = 0.0;

    printf("\n初始化完成,開始 3D 模擬...\n");

#ifdef _OPENMP
    double start_time_omp = omp_get_wtime();
    printf("Starting 3D simulation with OpenMP acceleration...\n");
#else
    struct timeval run_t_start, run_t_end;
    gettimeofday(&run_t_start, NULL);
    printf("Starting 3D simulation on Single Core...\n");
#endif

    while (t < t_end + dt_output) {

        if (t >= next_output_time) {

            printf("output_time: %.4f\n", t);

            char filename[256];
            snprintf(filename,
                     sizeof(filename),
                     "%s/output_%04d.csv",
                     output_folder,
                     output_step);

            write_csv(&sph, filename);

            printf("Step %d | Time: %.4f | dt: %.6f | Output: %s\n",
                   step,
                   t,
                   sph.dt,
                   filename);

            if (time_log != NULL) {
                fprintf(time_log, "%d,%.10f\n", output_step, t);
                fflush(time_log);
            }

            output_step++;
            next_output_time += dt_output;

            if (t >= t_end) {
                break;
            }
        }

        double dt = step_leapfrog_kdk_xreflective_yzperiodic_3d(
            &sph,
            compute_timestep_signal_velocity_3d,
            compute_force_xreflective_yperiodic_zperiodic_3d);

        t += dt;

        printf("sod3d t: %.10f, dt: %.10f\n", t, dt);

        step++;
    }

    double elapsed_time = 0.0;

#ifdef _OPENMP
    elapsed_time = omp_get_wtime() - start_time_omp;
#else
    gettimeofday(&run_t_end, NULL);
    elapsed_time =
        (run_t_end.tv_sec - run_t_start.tv_sec)
        + (run_t_end.tv_usec - run_t_start.tv_usec) / 1000000.0;
#endif

    printf("====================================\n");
    printf("3D Simulation finished!\n");
    printf("Total Execution Time: %.3f seconds\n", elapsed_time);
    printf("Total Outputfile numbers: %d\n", output_step);
    printf("Mass of particles: %f\n", mass);
    printf("====================================\n");

    if (time_log != NULL) {
        fclose(time_log);
    }

    free_sph_system(&sph);

    return 0;
}