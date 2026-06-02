#ifndef IO_H
#define IO_H
/*
 * Write SPH particle data to a CSV file.
 *
 * Output columns:
 * id, x, y, z, vx, vy, vz, ax, ay, az, m, rho, P, u, h
 */

void write_csv(SPHSystem * sph ,const char * filename);


#endif
