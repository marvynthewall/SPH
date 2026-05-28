# include "sph_system.h"
# include "force.h"

// this is for 1-layer loop
// fi, const = 2.0 for 2D, 3.0 for 3D
void compute_pressure_soundspeed_fi(SPHSystem2D *sph){
  double GG1 = sqrt(GAMMA * (GAMMA-1));
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].pressure = (GAMMA - 1) * sph->particles[i].density * sph->particles[i].u;
    sph->particles[i].cs = GG1 * sqrt(sph->particles[i].u);
    sph->particles[i].fi = 1.0 / (1 + sph->particles[i].h / (2.0 * sph->particles[i].density) * sph->particles[i].ddensity_dh);
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

  // initialize the ax ay in every particle
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
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
      double r_sqr = dx * dx + dy * dy;
      double r = sqrt(r_sqr);  // (optional?), sqrt is expensive

      double max_h = fmax(p_i->h, p_j->h);

      // Viscosity is active when r < max_h => r^2 < max_h^2
      if (r_sqr <= max_h * max_h) {

        // Pressure force
        // also a symmetry force
        // need to be < max_h ?
        // or only h?
        // F_P += ?
        
        


        // Viscosity force
        double f_v = 0.0;
        // as long as r < max_h, calculate it
        double dvx = p_i->vx - p_j->vx;
        double dvy = p_i->vy - p_j->vy;
        double r_dot_v = dx * dvx + dy * dvy;
        // Viscosity happens only when both particles are moving toward each other
        if (r_dot_v < 0.0){
          // Eq 31 mu_ij calculation
          double h_ij = (p_i->h + p_j->h) / 2.0;
          double mu_ij = h_ij * r_dot_v / (r_sqr + epsilon * (h_ij * h_ij));

          // Eq 30 PI_ij calculation
          double c_ij = (p_i->cs + p_j->cs) / 2.0;
          double rho_ij = (p_i->density + p_j->density) / 2.0;
          double PI_ij = (-alpha * c_ij * mu_ij + beta * mu_ij * mu_ij) / rho_ij;

          // Eq 29
          // du/dt(i)

          double W = 0.0, dWdr = 0.0, dWdh;
          cubic_spline_kernel_2d(r, p_i->h, &W, &dWdr, &dWdh);
          // Eq 27
          // W_ij avg

          // Eq 26
          // dv/dt(i)
          
          //F_V += ?
        }
      }

    }
  }
  return;
}

