# include "sph_system.h"


int check_particle_number(const SPHSystem2D *sph,
                                  int nx,
                                  int ny,
                                  const char *func_name)
{
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

/*
 * Initialize particles in a uniform 2D box.
 *
 * Particles are placed on a regular nx by ny grid.
 * This function only sets initial particle data.
 */
void init_uniform_box(SPHSystem2D *sph, int nx, int ny){
    int N_expected = check_particle_number(sph, nx, ny, "init_uniform_box");
    
    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;

    double dx = (xmax - xmin) / (nx-1);
    double dy = (ymax - ymin) / (ny-1);

    double mass = 1.0 / N_expected;
    double rho0 = 1.0;
    double P0   = 1.0;
    double u0   = 1.0;

    double h0 = 1.3 * dx;
    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
            int id = iy * nx + ix;
            /*
             * Place particles at cell centers.
             */
            sph->x[id] = xmin + (ix) * dx;
            sph->y[id] = ymin + (iy) * dy;

            /*
             * Initial velocity.
             */
            sph->vx[id] = 0.0;
            sph->vy[id] = 0.0;

            /*
             * Initial acceleration.
             */
            sph->ax[id] = 0.0;
            sph->ay[id] = 0.0;

            /*
             * Initial physical quantities.
             */
            sph->m[id]   = mass;
            sph->rho[id] = rho0;
            sph->P[id]   = P0;
            sph->u[id]   = u0;
            sph->h[id]   = h0;
        }
    }
}



void init_sod_2d(SPHSystem2D *sph, int nx, int ny){
    check_particle_number(sph, nx, ny, "init_sod_2d");
    
    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    double gamma = 1.4;

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

            sph->x[id] = x;
            sph->y[id] = y;

            sph->vx[id] = 0.0;
            sph->vy[id] = 0.0;

            sph->ax[id] = 0.0;
            sph->ay[id] = 0.0;

            if(x > 0.5){
                sph->rho[id] = rho_L;
                sph->P[id]   = P_L;
                sph->u[id]   = P_L / ((gamma - 1.0) * rho_L);
                sph->m[id]   = rho_L * dx * dy;
            }else{
                sph->rho[id] = rho_R;
                sph->P[id]   = P_R;
                sph->u[id]   = P_R / ((gamma - 1.0) * rho_R);
                sph->m[id]   = rho_R * dx * dy;
            }
            sph->h[id] = h0;
        }
    }
}