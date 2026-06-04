#ifndef RIEMANN_H
#define RIEMANN_H

#include "grid_system.h"

// HLLC Riemann solver
// W_L, W_R: Left and Right primitive variables at interface
// F: Computed flux (conservative) at interface
// dir: 0 for x-direction, 1 for y-direction
void hllc_flux(PrimVar *W_L, PrimVar *W_R, ConsVar *F, double gamma, int dir);

#endif
