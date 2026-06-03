#ifndef GRID_SYSTEM_H
#define GRID_SYSTEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double rho;      // Density
    double px, py;   // Momentum in x, y
    double E;        // Total energy per unit volume
} ConsVar;

typedef struct {
    double rho;      // Density
    double vx, vy;   // Velocity in x, y
    double P;        // Pressure
} PrimVar;

typedef struct {
    int id;          // Cell ID (useful for CSV output matching)
    double x, y;     // Cell center coordinates
    
    ConsVar U;       // Conservative variables (rho, rho v, E)
    PrimVar W;       // Primitive variables (rho, v, P)
    
    ConsVar U_new;   // Next step conservative variables
} Cell;

typedef struct {
    int nx, ny;          // Number of active cells in x, y
    int ng;              // Number of ghost cells (e.g., 2 for MUSCL-Hancock)
    int nx_tot, ny_tot;  // Total cells including ghost cells (nx + 2*ng)
    
    double dx, dy;       // Cell sizes
    double xmin, xmax;
    double ymin, ymax;
    
    double gamma;        // Specific heat ratio
    double cfl;          // Courant number
    
    double t, t_end;     // Time control
    double dt;           // Current time step
    
    Cell *cells;         // 1D array for 2D grid: index = i + j*nx_tot
} GridSystem;

// Macro for 2D to 1D index mapping
#define IDX(i, j, nx_tot) ((i) + (j) * (nx_tot))

// Helper Functions
void allocate_grid(GridSystem *grid, int nx, int ny, int ng, double xmin, double xmax, double ymin, double ymax, double gamma, double cfl, double t_end);
void free_grid(GridSystem *grid);

void prim_to_cons(PrimVar *W, ConsVar *U, double gamma);
void cons_to_prim(ConsVar *U, PrimVar *W, double gamma);
double get_soundspeed(double gamma, double P, double rho);

#endif
