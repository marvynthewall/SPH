#include "../include/grid_system.h"

void prim_to_cons(PrimVar *W, ConsVar *U, double gamma) {
    U->rho = W->rho;
    U->px = W->rho * W->vx;
    U->py = W->rho * W->vy;
    U->E = W->P / (gamma - 1.0) + 0.5 * W->rho * (W->vx * W->vx + W->vy * W->vy);
}

void cons_to_prim(ConsVar *U, PrimVar *W, double gamma) {
    W->rho = U->rho;
    W->vx = U->px / U->rho;
    W->vy = U->py / U->rho;
    W->P = (gamma - 1.0) * (U->E - 0.5 * U->rho * (W->vx * W->vx + W->vy * W->vy));
}

double get_soundspeed(double gamma, double P, double rho) {
    return sqrt(gamma * P / rho);
}

void allocate_grid(GridSystem *grid, int nx, int ny, int ng, double xmin, double xmax, double ymin, double ymax, double gamma, double cfl, double t_end) {
    grid->nx = nx;
    grid->ny = ny;
    grid->ng = ng;
    grid->nx_tot = nx + 2 * ng;
    grid->ny_tot = ny + 2 * ng;
    
    grid->dx = (xmax - xmin) / nx;
    if (ny > 1) {
        grid->dy = (ymax - ymin) / ny;
    } else {
        grid->dy = 1.0; // dummy for 1D
    }
    
    grid->xmin = xmin;
    grid->xmax = xmax;
    grid->ymin = ymin;
    grid->ymax = ymax;
    
    grid->gamma = gamma;
    grid->cfl = cfl;
    
    grid->t = 0.0;
    grid->t_end = t_end;
    grid->dt = 0.0;
    
    grid->cells = (Cell*)malloc(grid->nx_tot * grid->ny_tot * sizeof(Cell));
    
    for (int j = 0; j < grid->ny_tot; j++) {
        for (int i = 0; i < grid->nx_tot; i++) {
            int idx = IDX(i, j, grid->nx_tot);
            grid->cells[idx].id = idx;
            grid->cells[idx].x = xmin + (i - ng + 0.5) * grid->dx;
            if (ny > 1) {
                grid->cells[idx].y = ymin + (j - ng + 0.5) * grid->dy;
            } else {
                grid->cells[idx].y = 0.0;
            }
            memset(&grid->cells[idx].W, 0, sizeof(PrimVar));
            memset(&grid->cells[idx].U, 0, sizeof(ConsVar));
            memset(&grid->cells[idx].U_new, 0, sizeof(ConsVar));
        }
    }
}

void free_grid(GridSystem *grid) {
    if (grid->cells) {
        free(grid->cells);
    }
}
