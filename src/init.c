#include "sph_system.h"

int check_particle_number(const SPHSystem2D *sph,
                          int nx,
                          int ny,
                          const char *func_name)
{
    if (sph == NULL || sph->particles == NULL) {
        fprintf(stderr,
                "Error in %s: SPH system is not allocated.\n",
                func_name);
        exit(EXIT_FAILURE);
    }

    int N_expected = nx * ny;

    if (sph->N != N_expected) {
        fprintf(stderr,
                "Error in %s: sph->N = %d, but nx * ny = %d\n",
                func_name,
                sph->N,
                N_expected);
        exit(EXIT_FAILURE);
    }

    return N_expected;
}

void init_uniform_box(SPHSystem2D *sph, int nx, int ny)
{
    int N_expected = check_particle_number(sph, nx, ny, "init_uniform_box");

    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    double mass = 1.0 / N_expected;
    double rho0 = 1.0;
    double P0   = 1.0;
    double u0   = 1.0;
    double h0   = 1.3 * dx;

    double gamma = 1.4;

    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {

            int id = iy * nx + ix;

            sph->particles[id].id = id;

            sph->particles[id].x = xmin + ix * dx;
            sph->particles[id].y = ymin + iy * dy;

            sph->particles[id].vx = 0.0;
            sph->particles[id].vy = 0.0;

            sph->particles[id].ax = 0.0;
            sph->particles[id].ay = 0.0;

            sph->particles[id].mass     = mass;
            sph->particles[id].density  = rho0;
            sph->particles[id].pressure = P0;
            sph->particles[id].u        = u0;
            sph->particles[id].h        = h0;

            sph->particles[id].cs = sqrt(gamma * P0 / rho0);
        }
    }
}

void init_sod_2d(SPHSystem2D *sph, int nx, int ny)
{
    check_particle_number(sph, nx, ny, "init_sod_2d");

    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    double mass = 1.0 / N_expected;
    double rho0 = 1.0;
    double P0   = 1.0;
    double u0   = 1.0;
    double h0   = 1.3 * dx;

    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {

            int id = iy * nx + ix;

            sph->particles[id].id = id;

            sph->particles[id].x = xmin + ix * dx;
            sph->particles[id].y = ymin + iy * dy;

            sph->particles[id].vx = 0.0;
            sph->particles[id].vy = 0.0;

            sph->particles[id].ax = 0.0;
            sph->particles[id].ay = 0.0;

            sph->particles[id].mass     = mass;
            sph->particles[id].density  = rho0;
            sph->particles[id].pressure = P0;
            sph->particles[id].u        = u0;
            sph->particles[id].h        = h0;

            sph->particles[id].cs = 0.0;
        }
    }
}

void init_sod_2d(SPHSystem2D *sph, int nx, int ny)
{
    check_particle_number(sph, nx, ny, "init_sod_2d");

    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    double gamma = 5/3;

    double rho_L = 1.0;
    double P_L   = 1.0;

    double rho_R = 0.125;
    double P_R   = 0.1;

    double h0 = 1.3 * dx;

    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {

            int id = iy * nx + ix;

            double x = xmin + ix * dx;
            double y = ymin + iy * dy;

            sph->particles[id].id = id;

            sph->particles[id].x = x;
            sph->particles[id].y = y;

            sph->particles[id].vx = 0.0;
            sph->particles[id].vy = 0.0;

            sph->particles[id].ax = 0.0;
            sph->particles[id].ay = 0.0;

            if (x > 0.5) {
                sph->particles[id].density  = rho_L;
                sph->particles[id].pressure = P_L;
                sph->particles[id].u        = P_L / ((gamma - 1.0) * rho_L);
                sph->particles[id].mass     = rho_L * dx * dy;
            } else {
                sph->particles[id].density  = rho_R;
                sph->particles[id].pressure = P_R;
                sph->particles[id].u        = P_R / ((gamma - 1.0) * rho_R);
                sph->particles[id].mass     = rho_R * dx * dy;
            }

            sph->particles[id].h = h0;

            sph->particles[id].cs =
                sqrt(gamma *
                     sph->particles[id].pressure /
                     sph->particles[id].density);
        }
    }
}