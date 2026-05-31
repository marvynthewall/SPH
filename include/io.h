#ifndef IO_H
#define IO_H
/*
 * Write SPH particle data to a CSV file.
 *
 * Output columns:
 * id, x, y, vx, vy, ax, ay, m, rho, P, u, h
 */

void write_csv(SPHSystem * sph ,const char * filename);


#endif
