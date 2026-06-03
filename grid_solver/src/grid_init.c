#include "../include/grid_init.h"
#include <math.h>

void init_sod_1d(GridSystem *grid) {
  double x_mid = (grid->xmax + grid->xmin) / 2.0;
  for (int j = 0; j < grid->ny_tot; j++) {
    for (int i = 0; i < grid->nx_tot; i++) {
      int idx = IDX(i, j, grid->nx_tot);
      Cell *c = &grid->cells[idx];

      if (c->x <= x_mid) {
        c->W.rho = 1.0;
        c->W.P = 1.0;
      } else {
        c->W.rho = 0.125;
        c->W.P = 0.1;
      }
      c->W.vx = 0.0;
      c->W.vy = 0.0;

      prim_to_cons(&c->W, &c->U, grid->gamma);
    }
  }
}

void init_kh_2d(GridSystem *grid) {
  double rho_1 = 1.0;
  double rho_2 = 2.0;
  double v_1 = -0.5;
  double v_2 = 0.5;
  double P_0 = 2.5;
  double sigma = 0.025;

  for (int j = 0; j < grid->ny_tot; j++) {
    for (int i = 0; i < grid->nx_tot; i++) {
      int idx = IDX(i, j, grid->nx_tot);
      Cell *c = &grid->cells[idx];

      double px = c->x;
      double py = c->y;

      double f = 1.0 / ((1.0 + exp(-2.0 * (py - 0.25) / sigma)) *
                        (1.0 + exp(2.0 * (py - 0.75) / sigma)));

      c->W.rho = rho_1 + (rho_2 - rho_1) * f;
      c->W.vx = v_1 + (v_2 - v_1) * f;
      c->W.vy = 0.01 * sin(4.0 * M_PI * px);
      c->W.P = P_0;

      prim_to_cons(&c->W, &c->U, grid->gamma);
    }
  }
}
