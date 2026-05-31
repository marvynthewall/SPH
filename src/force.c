#include "force.h"
#include "sph_system.h"

// this is for 1-layer loop
// factor, const = 2.0 for 2D, 3.0 for 3D
// Eq(23)
void compute_pressure_soundspeed_factor(SPHSystem *sph) {
  double GG1 = sqrt(GAMMA * (GAMMA - 1));
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].pressure =
        (GAMMA - 1) * sph->particles[i].rho * sph->particles[i].u;
    sph->particles[i].cs = GG1 * sqrt(sph->particles[i].u);
    sph->particles[i].factor =
        1.0 / (1 + sph->particles[i].h / (2.0 * sph->particles[i].rho) *
                       sph->particles[i].drho_dh);
  }
  return;
}


// make sure that this is always inline
// inline: save time for jumping around different function part. 
// static: only inside this file, while compiling, this is better for optimization
// __attribute__((always_inline))
// when compiling, compiler can discard this after compiling this part, and optimize it directly. 
// inline is only a "suggestion", especially in O3, sometimes it changes idea, since this is a long function
// however in our simulation, this need to be inline!
__attribute__((always_inline)) static inline void compute_pairwise_physics(Particle *p_i, Particle *p_j, SPHSystem *sph) {
    double dx = p_i->x - p_j->x;
    double dy = p_i->y - p_j->y;
    double r = sqrt(dx * dx + dy * dy);
    if (r < 1e-12) return;  // avoid error

    double max_h = fmax(p_i->h, p_j->h);
    if (r > max_h) return;  // Viscosity & Force only active within max_h

    double dvx = p_i->vx - p_j->vx;
    double dvy = p_i->vy - p_j->vy;

    // get the W_ij(h_i) and W_ij(h_j)
    double W_i, dWdr_i, dWdh_i;
    cubic_spline_kernel_2d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
    double W_j, dWdr_j, dWdh_j;
    cubic_spline_kernel_2d(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

    // --- 1. Pressure Force ---
    // Pressure force, Eq(22)
    double term_i = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
    double term_j = p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;
    double scalar_force = p_j->mass * (term_i + term_j);

    double ax = -scalar_force * (dx/r);
    double ay = -scalar_force * (dy/r);

    p_i->ax += ax;
    p_i->ay += ay;

    // --- 2. Thermal Energy Evolution (Pressure Part) ---
    // update the specific thermal energy evolution, Eq(25)
    double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
    double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r); // 必須保留，avg_inner 會用到
    
    p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_j->mass * inner_product_v_dW_i;

    // --- 3. Artificial Viscosity ---
    double r_dot_v = dx * dvx + dy * dvy;
    double Vax=0.0, Vay=0.0;
    
    if (r_dot_v < 0.0) {
        // Eq 31 mu_ij calculation
        double h_ij = (p_i->h + p_j->h) / 2.0;
        double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
        // Eq 30 PI_ij calculation
        double c_ij = (p_i->cs + p_j->cs) / 2.0;
        double rho_ij = (p_i->rho + p_j->rho) / 2.0;
        double PI_ij = (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

        // update the specific thermal energy evolution
        // Eq(29)
        // inner product with the avg W is the avg of inner products
        double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
        p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

        // update the acceleration
        // Eq(26)
        double Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dx/r);
        double Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dy/r);
        p_i->ax += Vax;
        p_i->ay += Vay;
#ifndef OMP
        p_j->dudt += p_i->mass / 2.0 * PI_ij * avg_inner;
#endif
    }

    // --- 4. Reaction force only for single core
#ifndef OMP
    double mass_ratio = p_i->mass / p_j->mass;
    p_j->ax -= (ax + Vax) * mass_ratio;
    p_j->ay -= (ay + Vay) * mass_ratio;
    p_j->dudt += p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * p_i->mass * inner_product_v_dW_j;
#endif
}


void compute_force(SPHSystem *sph) {
#ifdef OMP
#pragma omp parallel for
#endif
  // initialize the ax,ay,dudt in every particle
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }
#ifdef OMP
#pragma omp parallel for
    // openmp, run every i to j, for parallelization
    for (int i = 0; i < sph->N; i++) {
        Particle *p_i = &sph->particles[i];
        for (int j = 0; j < sph->N; j++) {
            if (i == j) continue;
            Particle *p_j = &sph->particles[j];
            compute_pairwise_physics(p_i, p_j, sph);
        }
    }
#else
    // run only pairs (half the number of the computation)
    for (int i = 0; i < sph->N; i++) {
        Particle *p_i = &sph->particles[i];
        for (int j = i + 1; j < sph->N; j++) {
            Particle *p_j = &sph->particles[j];
            compute_pairwise_physics(p_i, p_j, sph);
        }
    }
#endif
}

// Computing force, symmetry 2-layers loop
void compute_force_xreflective_yperiodic(SPHSystem *sph) {
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];

    for (int j = i + 1; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j];

      double dy = p_i->y - p_j->y;
      if (dy > 0.5 * sph->box_size_y) dy -= sph->box_size_y;
      if (dy < -0.5 * sph->box_size_y) dy += sph->box_size_y;

      double dxs[3] = {p_i->x - p_j->x, p_i->x + p_j->x, p_i->x + p_j->x - 2.0 * sph->box_size_x};
      double dvxs[3] = {p_i->vx - p_j->vx, p_i->vx + p_j->vx, p_i->vx + p_j->vx};
      int signs_j[3] = {1, -1, -1}; // The reaction force on j has opposite sign for mirror interactions

      double dvy = p_i->vy - p_j->vy;
      double max_h = fmax(p_i->h, p_j->h);
      double mass_ratio = p_i->mass / p_j->mass;

      for (int k = 0; k < 3; k++) {
        double dx = dxs[k];
        double dvx = dvxs[k];
        int sign_j = signs_j[k];
        
        double r = sqrt(dx * dx + dy * dy);
        if (r < 1e-12 || r > max_h) continue;

        double W_i, dWdr_i, dWdh_i;
        cubic_spline_kernel_2d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
        double W_j, dWdr_j, dWdh_j;
        cubic_spline_kernel_2d(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

        double term_i = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
        double term_j = p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;

        double scalar_force = p_j->mass * (term_i + term_j);

        double ax = -scalar_force * (dx / r);
        double ay = -scalar_force * (dy / r);

        p_i->ax += ax;
        p_i->ay += ay;
        p_j->ax -= ax * mass_ratio * sign_j;
        p_j->ay -= ay * mass_ratio;

        double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
        double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r);
        p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_j->mass * inner_product_v_dW_i;
        p_j->dudt += p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * p_i->mass * inner_product_v_dW_j;

        double r_dot_v = dx * dvx + dy * dvy;
        if (r_dot_v < 0.0) {
          double h_ij = (p_i->h + p_j->h) / 2.0;
          double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
          double c_ij = (p_i->cs + p_j->cs) / 2.0;
          double rho_ij = (p_i->rho + p_j->rho) / 2.0;
          double PI_ij = (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

          double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
          p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;
          p_j->dudt += p_i->mass / 2.0 * PI_ij * avg_inner;

          double Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dx / r);
          double Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dy / r);
          p_i->ax += Vax;
          p_i->ay += Vay;
          p_j->ax -= Vax * mass_ratio * sign_j;
          p_j->ay -= Vay * mass_ratio;
        }
      }
    }
  }

  // Self-mirror interactions
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    double dxs[2] = {2.0 * p_i->x, 2.0 * (p_i->x - sph->box_size_x)};
    double dvxs[2] = {2.0 * p_i->vx, 2.0 * p_i->vx};

    for (int k = 0; k < 2; k++) {
      double dx = dxs[k];
      double dvx = dvxs[k];
      double r = fabs(dx);

      if (r < 1e-12 || r > p_i->h) continue;

      double W, dWdr, dWdh;
      cubic_spline_kernel_2d(r, p_i->h, &W, &dWdr, &dWdh);

      double term = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr;
      double scalar_force = p_i->mass * (2.0 * term);
      double ax = -scalar_force * (dx / r);

      p_i->ax += ax;

      double inner_product_v_dW = dvx * (dWdr * dx / r);
      p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_i->mass * inner_product_v_dW;

      double r_dot_v = dx * dvx;
      if (r_dot_v < 0.0) {
        double mu_ij = p_i->h * r_dot_v / (r * r + sph->epsilon * (p_i->h * p_i->h));
        double PI_ij = (-sph->alpha * p_i->cs * mu_ij + sph->beta * mu_ij * mu_ij) / p_i->rho;

        p_i->dudt += p_i->mass / 2.0 * PI_ij * inner_product_v_dW;
        double Vax = -p_i->mass * PI_ij * (dWdr * dx / r);
        p_i->ax += Vax;
      }
    }
  }
}

void compute_force_xperiodic_yperiodic(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {
      if (i == j) continue;
      Particle *p_j = &sph->particles[j];

      double dx = p_i->x - p_j->x;
      if (dx > 0.5 * sph->box_size_x) dx -= sph->box_size_x;
      else if (dx < -0.5 * sph->box_size_x) dx += sph->box_size_x;

      double dy = p_i->y - p_j->y;
      if (dy > 0.5 * sph->box_size_y) dy -= sph->box_size_y;
      else if (dy < -0.5 * sph->box_size_y) dy += sph->box_size_y;

      double r = sqrt(dx * dx + dy * dy);
      double max_h = fmax(p_i->h, p_j->h);
      
      if (r < 1e-12 || r > max_h) continue;

      double dvx = p_i->vx - p_j->vx;
      double dvy = p_i->vy - p_j->vy;

      double W_i, dWdr_i, dWdh_i;
      cubic_spline_kernel_2d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
      double W_j, dWdr_j, dWdh_j;
      cubic_spline_kernel_2d(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

      double term_i = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
      double term_j = p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;

      double scalar_force = p_j->mass * (term_i + term_j);

      double ax = -scalar_force * (dx / r);
      double ay = -scalar_force * (dy / r);

      p_i->ax += ax;
      p_i->ay += ay;

      double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
      double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r);
      p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_j->mass * inner_product_v_dW_i;

      double r_dot_v = dx * dvx + dy * dvy;
      if (r_dot_v < 0.0) {
        double h_ij = (p_i->h + p_j->h) / 2.0;
        double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
        double c_ij = (p_i->cs + p_j->cs) / 2.0;
        double rho_ij = (p_i->rho + p_j->rho) / 2.0;
        double PI_ij = (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

        double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
        p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

        double Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dx / r);
        double Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dy / r);
        p_i->ax += Vax;
        p_i->ay += Vay;
      }
    }
  }
}
