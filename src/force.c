# include "sph_system.h"
# include "force.h"

// this is for 1-layer loop
// factor, const = 2.0 for 2D, 3.0 for 3D
// Eq(23)
void compute_pressure_soundspeed_factor(SPHSystem2D *sph){
  double GG1 = sqrt(GAMMA * (GAMMA-1));
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].pressure = (GAMMA - 1) * sph->particles[i].density * sph->particles[i].u;
    sph->particles[i].cs = GG1 * sqrt(sph->particles[i].u);
    sph->particles[i].factor = 1.0 / (1 + sph->particles[i].h / (2.0 * sph->particles[i].density) * sph->particles[i].ddensity_dh);
  }
  return;
}

// Computing force, symmetry 2-layers loop
void compute_force(SPHSystem2D *sph){
  // for Viscositvy 
  // epsilon is usually set to 0.01, to prevent the too-close issue.
  // 0.5 <= alpha  <= 1.0
  // beta ~= 2 * alpha
  double epsilon = 0.01;
  double alpha = 1.0;
  double beta = 2.0 * alpha;

  // initialize the ax,ay,dudt in every particle
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

  // TODO spatial grid !!! can be used in every loop, at least also for the density

  // for every particle, find the surrounding and calculate
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i]; // indicate i-th particle

    // for every neighbor of the particle
    // because of symmetry, only calculate pairs
    for (int j = i+1; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j]; // indicate j-th particle

      // basic info between two particles 
      double dx = p_i->x - p_j->x;
      double dy = p_i->y - p_j->y;
      double r = sqrt(dx * dx + dy * dy);
      if (r < 1e-12) continue;  // avoid error

      double max_h = fmax(p_i->h, p_j->h);

      // Viscosity is active when r < max_h
      if (r <= max_h) {
        // some further basic info between two particles
        double dvx = p_i->vx - p_j->vx;
        double dvy = p_i->vy - p_j->vy;

        // get the W_ij(h_i) and W_ij(h_j)
        double W_i, dWdr_i, dWdh_i;
        cubic_spline_kernel_2d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
        double W_j, dWdr_j, dWdh_j;
        cubic_spline_kernel_2d(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

        // Pressure force, Eq(22)
        double term_i = p_i->factor * p_i->pressure / (p_i->density * p_i->density) * dWdr_i;
        double term_j = p_j->factor * p_j->pressure / (p_j->density * p_j->density) * dWdr_j;
        
        double scalar_force = p_j->mass * (term_i + term_j);
        
        double ax = -scalar_force * (dx/r);
        double ay = -scalar_force * (dy/r);
        
        // adding directly in the particle acceleration
        p_i->ax += ax;
        p_i->ay += ay;

        // the reaction force, negative direction, and consider the mass ratio
        double mass_ratio = p_i->mass / p_j->mass;
        p_j->ax -= ax * mass_ratio;
        p_j->ay -= ay * mass_ratio;

        // update the specific thermal energy evolution, Eq(25)
        // dv (dot) gradient of W(r, h_i) W(r, h_j)
        double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
        double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r);
        p_i->dudt += p_i->factor * p_i->pressure / (p_i->density * p_i->density) * p_j->mass * inner_product_v_dW_i;
        p_j->dudt += p_j->factor * p_j->pressure / (p_j->density * p_j->density) * p_i->mass * inner_product_v_dW_j;

        // Viscosity force
        // as long as r < max_h, calculate it
        double r_dot_v = dx * dvx + dy * dvy;
        // Viscosity happens only when both particles are moving toward each other
        if (r_dot_v < 0.0){
          // Eq 31 mu_ij calculation
          double h_ij = (p_i->h + p_j->h) / 2.0;
          double mu_ij = h_ij * r_dot_v / (r * r + epsilon * (h_ij * h_ij));

          // Eq 30 PI_ij calculation
          double c_ij = (p_i->cs + p_j->cs) / 2.0;
          double rho_ij = (p_i->density + p_j->density) / 2.0;
          double PI_ij = (-alpha * c_ij * mu_ij + beta * mu_ij * mu_ij) / rho_ij;

          // update the specific thermal energy evolution
          // Eq(29)
          // inner product with the avg W is the avg of inner products
          double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
          p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;
          p_j->dudt += p_i->mass / 2.0 * PI_ij * avg_inner;

          // update the acceleration
          // Eq(26)
          double Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dx/r);
          double Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dy/r);
          p_i->ax += Vax;
          p_i->ay += Vay;
          // reaction force
          p_j->ax -= Vax * mass_ratio;
          p_j->ay -= Vay * mass_ratio;
        }
      }
    }
  }
  return;
}

