#include "sph_system.h"

int check_particle_number(const SPHSystem *sph, int nx, int ny,
                          const char *func_name) {
  if (sph == NULL || sph->particles == NULL) {
    fprintf(stderr, "Error in %s: SPH system is not allocated.\n", func_name);
    exit(EXIT_FAILURE);
  }

  int N_expected = nx * ny;

  if (sph->N != N_expected) {
    fprintf(stderr, "Error in %s: sph->N = %d, but nx * ny = %d\n", func_name,
            sph->N, N_expected);
    exit(EXIT_FAILURE);
  }

  return N_expected;
}

void init_uniform_box(SPHSystem *sph, int nx, int ny) {
  int N_expected = check_particle_number(sph, nx, ny, "init_uniform_box");

  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;

  double dx = (xmax - xmin) / nx;
  double dy = (ymax - ymin) / ny;

  double mass = 1.0 / N_expected;
  double rho0 = 1.0;
  double P0 = 1.0;
  double u0 = 1.0;
  double h0 = 1.3 * dx;

  double gamma = 1.4;

  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {

      int id = iy * nx + ix;

      sph->particles[id].id = id;

      sph->particles[id].x = xmin + ((double)ix + 0.5) * dx;
      sph->particles[id].y = ymin + ((double)iy + 0.5) * dy;

      sph->particles[id].vx = 0.0;
      sph->particles[id].vy = 0.0;

      sph->particles[id].ax = 0.0;
      sph->particles[id].ay = 0.0;

      sph->particles[id].mass = mass;
      sph->particles[id].rho = rho0;
      sph->particles[id].pressure = P0;
      sph->particles[id].u = u0;
      sph->particles[id].h = h0;

      sph->particles[id].cs = sqrt(gamma * P0 / rho0);
    }
  }
}

void init_sod_2d(SPHSystem *sph, int nx, int ny) {
  check_particle_number(sph, nx, ny, "init_sod_2d");

  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;

  double dx = (xmax - xmin) / (nx - 1);
  double dy = (ymax - ymin) / (ny - 1);

  double gamma = 1.4;

  double rho_L = 1.0;
  double P_L = 1.0;

  double rho_R = 0.125;
  double P_R = 0.1;

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
        sph->particles[id].rho = rho_L;
        sph->particles[id].pressure = P_L;
        sph->particles[id].u = P_L / ((gamma - 1.0) * rho_L);
        sph->particles[id].mass = rho_L * dx * dy;
      } else {
        sph->particles[id].rho = rho_R;
        sph->particles[id].pressure = P_R;
        sph->particles[id].u = P_R / ((gamma - 1.0) * rho_R);
        sph->particles[id].mass = rho_R * dx * dy;
      }

      sph->particles[id].h = h0;

      sph->particles[id].cs =
          sqrt(gamma * sph->particles[id].pressure / sph->particles[id].rho);
    }
  }
}

void init_sph_parameter(SPHSystem *sph) {
  // Viscosity
  // epsilon is usually set to 0.01, to prevent the too-close issue.
  // 0.5 <= alpha  <= 1.0
  // beta ~= 2 * alpha
  sph->epsilon = 0.01;
  sph->alpha = 1.0;
  sph->beta = 2.0 * sph->alpha;
}

void calculate_position(int N, double *pos_x, double *pos_y, double x_min,
                        double x_max, double y_min, double y_max, char tail) {
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

  while (nx * ny < N) {
    if (ny < (int)floor(r * (double)nx))
      ny += 1;
    else
      nx += 1;
  }

  int n_tail = N - (nx - 1) * ny;
  // distribute the position
  int n = 0;
  double dx = Lx / (double)(nx);
  double dy = Ly / (double)(ny);
  for (int i = 0; i < nx; i++) {
    if ((tail == 'l' && i == 0) || (tail == 'r' && i == nx - 1)) {
      double tdy = Ly / (double)n_tail;
      for (int j = 0; j < n_tail; j++) {
        pos_x[n] = x_min + ((double)i + 0.5) * dx;
        pos_y[n] = y_min + ((double)j + 0.5) * tdy;
        n++;
      }
    } else {
      for (int j = 0; j < ny; j++) {
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
void init_sod_2d_2(SPHSystem *sph, double x, double y, double mass) {
  // check_particle_number(sph, nx, ny, "init_sod_2d");
  double gamma = 1.4;

  double rho_L = 1.0;
  double P_L = 1.0;

  double rho_R = 0.125;
  double P_R = 0.1;

  // mass per particle
  // double mass = 0.001;

  double eta = 1.3;
  double h_L = sqrt(eta * eta * mass / rho_L);
  double h_R = sqrt(eta * eta * mass / rho_R);

  double x_mid = x / 2.0;

  // left 0 ~ x_mid, right: x_mid ~ x
  double total_mass_L = rho_L * x_mid * y;
  double total_mass_R = rho_R * x_mid * y;
  double total_mass = total_mass_L + total_mass_R;
  int N = (int)round(total_mass / mass);
  int N_L = (int)round(total_mass_L / mass);
  int N_R = N - N_L;

  allocate_sph_system(sph, N);
  sph->gamma = gamma;

  sph->box_size_x = x;
  sph->box_size_y = y;

  // ID for left: 0 ~ N_L-1
  // ID for right: N_L ~ N-1

  double pos_x_L[N_L], pos_y_L[N_L];
  double pos_x_R[N_R], pos_y_R[N_R];

  calculate_position(N_L, pos_x_L, pos_y_L, 0, x_mid, 0, y, 'l');
  calculate_position(N_R, pos_x_R, pos_y_R, x_mid, x, 0, y, 'r');

  for (int i = 0; i < N_L; i++) {
    sph->particles[i].id = i;

    sph->particles[i].x = pos_x_L[i];
    sph->particles[i].y = pos_y_L[i];
    sph->particles[i].rho = rho_L;
    sph->particles[i].pressure = P_L;
    sph->particles[i].u = P_L / ((gamma - 1.0) * rho_L);
    sph->particles[i].h = h_L;
    sph->particles[i].cs = sqrt(gamma * P_L / rho_L);
  }
  for (int i = N_L; i < N; i++) {
    sph->particles[i].id = i;

    sph->particles[i].x = pos_x_R[i - N_L];
    sph->particles[i].y = pos_y_R[i - N_L];
    sph->particles[i].rho = rho_R;
    sph->particles[i].pressure = P_R;
    sph->particles[i].u = P_R / ((gamma - 1.0) * rho_R);
    sph->particles[i].h = h_R;
    sph->particles[i].cs = sqrt(gamma * P_R / rho_R);
  }
}

void init_sod_2d_3(SPHSystem *sph, double x_max, double y_max,
                   double target_mass) {
  double gamma = 1.4;

  // Macroscopic physical quantities
  double rho_L = 1.0, P_L = 1.0;
  double rho_R = 0.125, P_R = 0.1;
  double eta = 1.3;

  double x_mid = x_max / 2.0;

  // 1. Derive the perfect square spacing (dx = dy) from target mass and density
  double dx_L = sqrt(target_mass / rho_L);
  double dx_R = sqrt(target_mass / rho_R);

  // 2. Calculate the number of particles each region can hold
  // (floor the value to ensure particles stay within the physical boundaries)
  int nx_L = (int)(x_mid / dx_L);
  int ny_L = (int)(y_max / dx_L);

  int nx_R = (int)((x_max - x_mid) / dx_R);
  int ny_R = (int)(y_max / dx_R);

  // 3. Let the geometric space determine the actual total number of particles
  int N_L = nx_L * ny_L;
  int N_R = nx_R * ny_R;
  int N = N_L + N_R;

  printf("Number of Particles: %d\n", N);

  // Allocate memory for the SPH system
  allocate_sph_system(sph, N);

  sph->gamma = gamma;
  sph->box_size_x = x_max;
  sph->box_size_y = y_max;

  // Pre-calculate smoothing lengths for both regions
  double h_L = eta * dx_L;
  double h_R = eta * dx_R;

  int p_idx = 0; // Particle array index

  // 4. Distribute the high-pressure fluid in the left region
  // (Center the grid by adding an offset, leaving half a dx as the distance to
  // the walls)
  double offset_x_L = (x_mid - (nx_L * dx_L)) / 2.0 + (dx_L / 2.0);
  double offset_y_L = (y_max - (ny_L * dx_L)) / 2.0 + (dx_L / 2.0);

  for (int i = 0; i < nx_L; i++) {
    for (int j = 0; j < ny_L; j++) {
      Particle *p = &sph->particles[p_idx++];
      p->id = p_idx;
      p->x = offset_x_L + i * dx_L;
      p->y = offset_y_L + j * dx_L;
      p->mass = target_mass;
      p->rho = rho_L;
      p->pressure = P_L;
      p->u = P_L / ((gamma - 1.0) * rho_L);
      p->h = h_L;
      p->cs = sqrt(gamma * P_L / rho_L);
      p->vx = 0.0;
      p->vy = 0.0;
    }
  }

  // 5. Distribute the low-pressure fluid in the right region
  // (X coordinates start from x_mid)
  double offset_x_R =
      x_mid + ((x_max - x_mid - (nx_R * dx_R)) / 2.0) + (dx_R / 2.0);
  double offset_y_R = (y_max - (ny_R * dx_R)) / 2.0 + (dx_R / 2.0);

  for (int i = 0; i < nx_R; i++) {
    for (int j = 0; j < ny_R; j++) {
      Particle *p = &sph->particles[p_idx++];
      p->id = p_idx;
      p->x = offset_x_R + i * dx_R;
      p->y = offset_y_R + j * dx_R;
      p->mass = target_mass;
      p->rho = rho_R;
      p->pressure = P_R;
      p->u = P_R / ((gamma - 1.0) * rho_R);
      p->h = h_R;
      p->cs = sqrt(gamma * P_R / rho_R);
      p->vx = 0.0;
      p->vy = 0.0;
    }
  }
}

void init_sod_3d_3(SPHSystem *sph, double x_max, double y_max, double z_max,
                   double target_mass) {
  double gamma = 1.4;

  // Macroscopic physical quantities
  double rho_L = 1.0, P_L = 1.0;
  double rho_R = 0.125, P_R = 0.1;
  double eta = 1.3;

  double x_mid = x_max / 2.0;

  // 1. Derive cubic spacing from target mass and density
  double dx_L = cbrt(target_mass / rho_L);
  double dx_R = cbrt(target_mass / rho_R);

  // 2. Calculate number of particles in each region
  int nx_L = (int)(x_mid / dx_L);
  int ny_L = (int)(y_max / dx_L);
  int nz_L = (int)(z_max / dx_L);

  int nx_R = (int)((x_max - x_mid) / dx_R);
  int ny_R = (int)(y_max / dx_R);
  int nz_R = (int)(z_max / dx_R);

  // 3. Total particle number
  int N_L = nx_L * ny_L * nz_L;
  int N_R = nx_R * ny_R * nz_R;
  int N = N_L + N_R;

  printf("Number of Particles: %d\n", N);

  allocate_sph_system(sph, N);

  sph->dim = 3;
  sph->gamma = gamma;
  sph->box_size_x = x_max;
  sph->box_size_y = y_max;
  sph->box_size_z = z_max;

  // 4. Smoothing lengths
  double h_L = eta * dx_L;
  double h_R = eta * dx_R;

  int p_idx = 0;

  // 5. Left region: high density / high pressure
  double offset_x_L = (x_mid - (nx_L * dx_L)) / 2.0 + (dx_L / 2.0);
  double offset_y_L = (y_max - (ny_L * dx_L)) / 2.0 + (dx_L / 2.0);
  double offset_z_L = (z_max - (nz_L * dx_L)) / 2.0 + (dx_L / 2.0);

  for (int i = 0; i < nx_L; i++) {
    for (int j = 0; j < ny_L; j++) {
      for (int k = 0; k < nz_L; k++) {

        Particle *p = &sph->particles[p_idx++];

        p->id = p_idx;

        p->x = offset_x_L + i * dx_L;
        p->y = offset_y_L + j * dx_L;
        p->z = offset_z_L + k * dx_L;

        p->mass = target_mass;

        p->rho = rho_L;
        p->pressure = P_L;
        p->u = P_L / ((gamma - 1.0) * rho_L);

        p->h = h_L;
        p->cs = sqrt(gamma * P_L / rho_L);

        p->vx = 0.0;
        p->vy = 0.0;
        p->vz = 0.0;

        p->ax = 0.0;
        p->ay = 0.0;
        p->az = 0.0;
        p->drho_dh = 0.0;
        p->factor = 1.0;
      }
    }
  }

  // 6. Right region: low density / low pressure
  double offset_x_R =
      x_mid + ((x_max - x_mid - (nx_R * dx_R)) / 2.0) + (dx_R / 2.0);
  double offset_y_R = (y_max - (ny_R * dx_R)) / 2.0 + (dx_R / 2.0);
  double offset_z_R = (z_max - (nz_R * dx_R)) / 2.0 + (dx_R / 2.0);

  for (int i = 0; i < nx_R; i++) {
    for (int j = 0; j < ny_R; j++) {
      for (int k = 0; k < nz_R; k++) {

        Particle *p = &sph->particles[p_idx++];

        p->id = p_idx;

        p->x = offset_x_R + i * dx_R;
        p->y = offset_y_R + j * dx_R;
        p->z = offset_z_R + k * dx_R;

        p->mass = target_mass;

        p->rho = rho_R;
        p->pressure = P_R;
        p->u = P_R / ((gamma - 1.0) * rho_R);

        p->h = h_R;
        p->cs = sqrt(gamma * P_R / rho_R);

        p->vx = 0.0;
        p->vy = 0.0;
        p->vz = 0.0;

        p->ax = 0.0;
        p->ay = 0.0;
        p->az = 0.0;
      }
    }
  }
}

void init_KH(SPHSystem *sph, int nx, int ny) {
  // Use constants from section 5.5 Contact Discontinuities and Fluid
  // Instabilities
  double gamma = 5.0 / 3.0;
  double rho_1 = 1.0;
  double rho_2 = 2.0;
  double v_1 = -0.5;
  double v_2 = 0.5;
  double P_0 = 2.5;

  double sigma = 0.025;
  double eta = 1.8; // smoothing factor

  double x_max = 1.0;
  double y_max = 1.0;

  double dx = x_max / nx;
  double dy = y_max / ny;

  int N = nx * ny;
  allocate_sph_system(sph, N);

  sph->box_size_x = x_max;
  sph->box_size_y = y_max;

  double h_uniform = eta * sqrt(dx * dy);

  int p_idx = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      Particle *p = &sph->particles[p_idx++];
      p->id = p_idx;

      double px = (i + 0.5) * dx;
      double py = (j + 0.5) * dy;

      // Calculate smoothed profiles (eq. 71)
      double f = 1.0 / ((1.0 + exp(-2.0 * (py - 0.25) / sigma)) *
                        (1.0 + exp(2.0 * (py - 0.75) / sigma)));

      double rho = rho_1 + (rho_2 - rho_1) * f;
      double vx = v_1 + (v_2 - v_1) * f;
      // trigger instability by y-direction velocity perturbation
      double vy = 0.01 * sin(4.0 * M_PI * px);

      double mass = rho * dx * dy;

      p->x = px;
      p->y = py;
      p->vx = vx;
      p->vy = vy;
      p->mass = mass;
      p->rho = rho;
      p->pressure = P_0;
      p->u = P_0 / ((gamma - 1.0) * rho);
      p->h = h_uniform;
      p->cs = sqrt(gamma * P_0 / rho);
    }
  }
}

void init_sod_1d(SPHSystem *sph, double x_max, double target_mass) {
  double gamma = 1.4;

  double rho_L = 1.0, P_L = 1.0;
  double rho_R = 0.125, P_R = 0.1;
  double eta = 1.3;

  double x_mid = x_max / 2.0;

  double dx_L = target_mass / rho_L;
  double dx_R = target_mass / rho_R;

  int nx_L = (int)(x_mid / dx_L);
  int nx_R = (int)((x_max - x_mid) / dx_R);

  int N = nx_L + nx_R;
  printf("Number of Particles (1D): %d\n", N);

  allocate_sph_system(sph, N);

  sph->gamma = gamma;
  sph->box_size_x = x_max;
  sph->box_size_y = 1.0; // Dummy
  sph->box_size_z = 1.0; // Dummy

  double h_L = eta * dx_L;
  double h_R = eta * dx_R;

  int p_idx = 0;

  double offset_x_L = (x_mid - (nx_L * dx_L)) / 2.0 + (dx_L / 2.0);

  for (int i = 0; i < nx_L; i++) {
    Particle *p = &sph->particles[p_idx++];
    p->id = p_idx;
    p->x = offset_x_L + i * dx_L;
    p->y = 0.0;
    p->z = 0.0;
    p->mass = target_mass;
    p->rho = rho_L;
    p->pressure = P_L;
    p->u = P_L / ((gamma - 1.0) * rho_L);
    p->h = h_L;
    p->cs = sqrt(gamma * P_L / rho_L);
    p->vx = 0.0;
    p->vy = 0.0;
    p->vz = 0.0;
    p->ax = 0.0;
    p->ay = 0.0;
    p->az = 0.0;
    p->drho_dh = 0.0;
    p->factor = 1.0;
  }

  double offset_x_R =
      x_mid + ((x_max - x_mid - (nx_R * dx_R)) / 2.0) + (dx_R / 2.0);

  for (int i = 0; i < nx_R; i++) {
    Particle *p = &sph->particles[p_idx++];
    p->id = p_idx;
    p->x = offset_x_R + i * dx_R;
    p->y = 0.0;
    p->z = 0.0;
    p->mass = target_mass;
    p->rho = rho_R;
    p->pressure = P_R;
    p->u = P_R / ((gamma - 1.0) * rho_R);
    p->h = h_R;
    p->cs = sqrt(gamma * P_R / rho_R);
    p->vx = 0.0;
    p->vy = 0.0;
    p->vz = 0.0;
    p->ax = 0.0;
    p->ay = 0.0;
    p->az = 0.0;
    p->drho_dh = 0.0;
    p->factor = 1.0;
  }
}
