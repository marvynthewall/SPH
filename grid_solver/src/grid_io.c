#include "../include/grid_io.h"

void write_csv(GridSystem *grid, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "id,x,y,z,vx,vy,vz,ax,ay,az,m,rho,P,u,h,cs\n");
    
    int out_id = 0;
    for (int j = grid->ng; j < grid->ny_tot - grid->ng; j++) {
        for (int i = grid->ng; i < grid->nx_tot - grid->ng; i++) {
            int idx = IDX(i, j, grid->nx_tot);
            Cell *c = &grid->cells[idx];
            
            double cell_vol = grid->dx * (grid->ny > 1 ? grid->dy : 1.0);
            double mass = c->W.rho * cell_vol;
            double u = c->W.P / ((grid->gamma - 1.0) * c->W.rho);
            double cs = get_soundspeed(grid->gamma, c->W.P, c->W.rho);
            double h = fmin(grid->dx, (grid->ny > 1 ? grid->dy : grid->dx)) * 1.3;
            
            fprintf(fp, "%d,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                    out_id++, c->x, c->y, 0.0,
                    c->W.vx, c->W.vy, 0.0,
                    0.0, 0.0, 0.0,
                    mass, c->W.rho, c->W.P, u, h, cs);
        }
    }
    
    fclose(fp);
}
