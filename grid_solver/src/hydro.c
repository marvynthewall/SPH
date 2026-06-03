#include "../include/hydro.h"
#include "../include/riemann.h"

// Minmod limiter
double minmod(double a, double b) {
  if (a * b <= 0)
    return 0.0;
  if (fabs(a) < fabs(b))
    return a;
  return b;
}

void calculate_dt(GridSystem *grid) {
  double max_speed_x = 0.0;
  double max_speed_y = 0.0;

  for (int j = grid->ng; j < grid->ny_tot - grid->ng; j++) {
    for (int i = grid->ng; i < grid->nx_tot - grid->ng; i++) {
      int idx = IDX(i, j, grid->nx_tot);
      Cell *c = &grid->cells[idx];

      double cs = get_soundspeed(grid->gamma, c->W.P, c->W.rho);
      double sx = fabs(c->W.vx) + cs;
      double sy = fabs(c->W.vy) + cs;

      if (sx > max_speed_x)
        max_speed_x = sx;
      if (sy > max_speed_y)
        max_speed_y = sy;
    }
  }

  double dt_x = grid->dx / max_speed_x;
  double dt_y = (grid->ny > 1) ? (grid->dy / max_speed_y) : 1e10;

  grid->dt = grid->cfl * fmin(dt_x, dt_y);
  if (grid->t + grid->dt > grid->t_end) {
    grid->dt = grid->t_end - grid->t;
  }
}

void apply_boundary_conditions(GridSystem *grid, int prob_type) {
  int ng = grid->ng;
  int nx_tot = grid->nx_tot;
  int ny_tot = grid->ny_tot;

  // 0 = Sod (Transmissive X), 1 = KH (Periodic X, Transmissive Y)

  if (prob_type == 0) {
    // Transmissive X
    for (int j = 0; j < ny_tot; j++) {
      for (int i = 0; i < ng; i++) {
        grid->cells[IDX(i, j, nx_tot)].W = grid->cells[IDX(ng, j, nx_tot)].W;
        grid->cells[IDX(i, j, nx_tot)].U = grid->cells[IDX(ng, j, nx_tot)].U;

        grid->cells[IDX(nx_tot - 1 - i, j, nx_tot)].W =
            grid->cells[IDX(nx_tot - 1 - ng, j, nx_tot)].W;
        grid->cells[IDX(nx_tot - 1 - i, j, nx_tot)].U =
            grid->cells[IDX(nx_tot - 1 - ng, j, nx_tot)].U;
      }
    }
  } else if (prob_type == 1) {
    // Periodic X
    for (int j = 0; j < ny_tot; j++) {
      for (int i = 0; i < ng; i++) {
        grid->cells[IDX(i, j, nx_tot)].W =
            grid->cells[IDX(nx_tot - 2 * ng + i, j, nx_tot)].W;
        grid->cells[IDX(i, j, nx_tot)].U =
            grid->cells[IDX(nx_tot - 2 * ng + i, j, nx_tot)].U;

        grid->cells[IDX(nx_tot - ng + i, j, nx_tot)].W =
            grid->cells[IDX(ng + i, j, nx_tot)].W;
        grid->cells[IDX(nx_tot - ng + i, j, nx_tot)].U =
            grid->cells[IDX(ng + i, j, nx_tot)].U;
      }
    }
    // Transmissive Y
    for (int i = 0; i < nx_tot; i++) {
      for (int j = 0; j < ng; j++) {
        grid->cells[IDX(i, j, nx_tot)].W = grid->cells[IDX(i, ng, nx_tot)].W;
        grid->cells[IDX(i, j, nx_tot)].U = grid->cells[IDX(i, ng, nx_tot)].U;
        // Reverse vy for reflective boundary (optional, usually transmissive is
        // ok for KH if box is large enough, let's use transmissive)

        grid->cells[IDX(i, ny_tot - 1 - j, nx_tot)].W =
            grid->cells[IDX(i, ny_tot - 1 - ng, nx_tot)].W;
        grid->cells[IDX(i, ny_tot - 1 - j, nx_tot)].U =
            grid->cells[IDX(i, ny_tot - 1 - ng, nx_tot)].U;
      }
    }
  }
}

// 1st order Godunov step for simplicity and robustness
void hydro_step(GridSystem *grid) {
  int nx_tot = grid->nx_tot;
  int ny_tot = grid->ny_tot;
  int ng = grid->ng;

  // Arrays for fluxes
  ConsVar *F_x = (ConsVar *)malloc(nx_tot * ny_tot * sizeof(ConsVar));
  ConsVar *F_y = (ConsVar *)malloc(nx_tot * ny_tot * sizeof(ConsVar));
  memset(F_x, 0, nx_tot * ny_tot * sizeof(ConsVar));
  memset(F_y, 0, nx_tot * ny_tot * sizeof(ConsVar));

  // Calculate X-fluxes at interfaces i+1/2
  for (int j = 0; j < ny_tot; j++) {
    for (int i = ng - 1; i < nx_tot - ng; i++) {
      int idx_L = IDX(i, j, nx_tot);
      int idx_R = IDX(i + 1, j, nx_tot);
      hllc_flux(&grid->cells[idx_L].W, &grid->cells[idx_R].W, &F_x[idx_L],
                grid->gamma, 0);
    }
  }

  // Calculate Y-fluxes at interfaces j+1/2
  if (grid->ny > 1) {
    for (int j = ng - 1; j < ny_tot - ng; j++) {
      for (int i = 0; i < nx_tot; i++) {
        int idx_L = IDX(i, j, nx_tot);
        int idx_R = IDX(i, j + 1, nx_tot);
        hllc_flux(&grid->cells[idx_L].W, &grid->cells[idx_R].W, &F_y[idx_L],
                  grid->gamma, 1);
      }
    }
  }

  // Update conservative variables
  double dtdx = grid->dt / grid->dx;
  double dtdy = (grid->ny > 1) ? (grid->dt / grid->dy) : 0.0;

  for (int j = ng; j < ny_tot - ng; j++) {
    for (int i = ng; i < nx_tot - ng; i++) {
      int idx = IDX(i, j, nx_tot);
      int idx_L = IDX(i - 1, j, nx_tot);
      int idx_B = IDX(i, j - 1, nx_tot);

      grid->cells[idx].U_new.rho =
          grid->cells[idx].U.rho - dtdx * (F_x[idx].rho - F_x[idx_L].rho);
      grid->cells[idx].U_new.px =
          grid->cells[idx].U.px - dtdx * (F_x[idx].px - F_x[idx_L].px);
      grid->cells[idx].U_new.py =
          grid->cells[idx].U.py - dtdx * (F_x[idx].py - F_x[idx_L].py);
      grid->cells[idx].U_new.E =
          grid->cells[idx].U.E - dtdx * (F_x[idx].E - F_x[idx_L].E);

      if (grid->ny > 1) {
        grid->cells[idx].U_new.rho -= dtdy * (F_y[idx].rho - F_y[idx_B].rho);
        grid->cells[idx].U_new.px -= dtdy * (F_y[idx].px - F_y[idx_B].px);
        grid->cells[idx].U_new.py -= dtdy * (F_y[idx].py - F_y[idx_B].py);
        grid->cells[idx].U_new.E -= dtdy * (F_y[idx].E - F_y[idx_B].E);
      }
    }
  }

  // Copy new values and update primitives
  for (int j = ng; j < ny_tot - ng; j++) {
    for (int i = ng; i < nx_tot - ng; i++) {
      int idx = IDX(i, j, nx_tot);
      grid->cells[idx].U = grid->cells[idx].U_new;
      cons_to_prim(&grid->cells[idx].U, &grid->cells[idx].W, grid->gamma);
    }
  }

  free(F_x);
  free(F_y);
}
