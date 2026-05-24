#ifndef INIT_H
#define INIT_H

int check_particle_number(const SPHSystem2D *sph,
                                  int nx,
                                  int ny,
                                  const char *func_name);

void init_uniform_box(SPHSystem2D *sph, int nx, int ny);

void init_sod_2d(SPHSystem2D *sph, int nx, int ny);

#endif