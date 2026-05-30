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
            sph->particles[id].rho      = rho0;
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

            sph->particles[id].id = id;

            sph->particles[id].x = x;
            sph->particles[id].y = y;

            sph->particles[id].vx = 0.0;
            sph->particles[id].vy = 0.0;

            sph->particles[id].ax = 0.0;
            sph->particles[id].ay = 0.0;

            if (x > 0.5) {
                sph->particles[id].rho  = rho_L;
                sph->particles[id].pressure = P_L;
                sph->particles[id].u        = P_L / ((gamma - 1.0) * rho_L);
                sph->particles[id].mass     = rho_L * dx * dy;
            } else {
                sph->particles[id].rho  = rho_R;
                sph->particles[id].pressure = P_R;
                sph->particles[id].u        = P_R / ((gamma - 1.0) * rho_R);
                sph->particles[id].mass     = rho_R * dx * dy;
            }

            sph->particles[id].h = h0;

            sph->particles[id].cs =
                sqrt(gamma *
                     sph->particles[id].pressure /
                     sph->particles[id].rho);
        }
    }
}

void init_sph_parameter(SPHSystem2D *sph){
    // Viscosity
    // epsilon is usually set to 0.01, to prevent the too-close issue.
    // 0.5 <= alpha  <= 1.0
    // beta ~= 2 * alpha
    sph->epsilon = 0.01;
    sph->alpha   = 1.0;
    sph->beta    = 2.0 * sph->alpha;
}

void calculate_position(int N, double* pos_x, double* pos_y, double x_min, double x_max, double y_min, double y_max, char tail){
    // tail = 'l' or 'r'
    // size
    double Lx = x_max - x_min;
    double Ly = y_max - y_min;
    // ratio
    double r = Ly / Lx;
    
    // in x direction, need around nx particles;
    // find nx
    int nx = (int)floor(sqrt((double)N / r));
    int ny = (int)floor(r * (double)nx);
    // now nx * ny < N
    
    while (nx * ny < N){
        if (ny < (int)floor(r * (double)nx))
            ny += 1;
        else
            nx += 1;
    }

    int n_tail = N - (nx-1) * ny;
    // distribute the position
    int n = 0;
    double dx = Lx / (double)(nx);
    double dy = Ly / (double)(ny);
    for(int i = 0 ; i < nx; i++){
        if((tail == 'l' && i == 0) || (tail == 'r' && i == nx-1)){
            double tdy = Ly / (double)n_tail;
            for(int j = 0; j < n_tail ; j++){
                pos_x[n] = x_min + ((double)i + 0.5) * dx;
                pos_y[n] = y_min + ((double)j + 0.5) * tdy;
                n++;
            }
        }
        else{
            for(int j = 0 ; j < ny; j++){
                pos_x[n] = x_min + ((double)i + 0.5) * dx;
                pos_y[n] = y_min + ((double)j + 0.5) * dy;
                n++;
            }
        }
    }
    // check error
    if (n != N)
        printf("calculate position n != N\n");
    return;
}

// x, y is the physical size of the boxes
void init_sod_2d_2(SPHSystem2D *sph, double x, double y, double mass)
{
    // check_particle_number(sph, nx, ny, "init_sod_2d");
    double gamma = 1.4;

    double rho_L = 1.0;
    double P_L   = 1.0;

    double rho_R = 0.125;
    double P_R   = 0.1;


    // mass per particle
    // double mass = 0.001;
    
    double eta = 1.3;
    double h_L = sqrt(eta*eta*mass/rho_L);
    double h_R = sqrt(eta*eta*mass/rho_R);

    double x_mid = x / 2.0;

    // left 0 ~ x_mid, right: x_mid ~ x
    double total_mass_L = rho_L * x_mid * y;
    double total_mass_R = rho_R * x_mid * y;
    double total_mass = total_mass_L + total_mass_R;
    int N = (int)round(total_mass / mass);
    int N_L = (int)round(total_mass_L / mass);
    int N_R = N - N_L;
    
    allocate_sph_system(sph, N);
    
    // ID for left: 0 ~ N_L-1
    // ID for right: N_L ~ N-1

    double pos_x_L[N_L], pos_y_L[N_L];
    double pos_x_R[N_R], pos_y_R[N_R];

    calculate_position(N_L, pos_x_L, pos_y_L, 0, x_mid, 0, y, 'l');
    calculate_position(N_R, pos_x_R, pos_y_R, x_mid, x, 0, y, 'r');

    for (int i = 0 ; i < N_L ; i++){
            sph->particles[i].id = i;

            sph->particles[i].x        = pos_x_L[i];
            sph->particles[i].y        = pos_y_L[i];
            sph->particles[i].rho      = rho_L;
            sph->particles[i].pressure = P_L;
            sph->particles[i].u        = P_L / ((gamma - 1.0) * rho_L);
            sph->particles[i].h = h_L;
            sph->particles[i].cs =sqrt(gamma * P_L / rho_L);
    }
    for (int i = N_L ; i < N ; i++){
            sph->particles[i].id = i;

            sph->particles[i].x        = pos_x_R[i-N_L];
            sph->particles[i].y        = pos_y_R[i-N_L];
            sph->particles[i].rho      = rho_R;
            sph->particles[i].pressure = P_R;
            sph->particles[i].u        = P_R / ((gamma - 1.0) * rho_R);
            sph->particles[i].h = h_R;
            sph->particles[i].cs =sqrt(gamma * P_R / rho_R);
    }
}
