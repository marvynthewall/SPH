#ifndef INIT_H
#define INIT_H

int check_particle_number(const SPHSystem *sph,
                                  int nx,
                                  int ny,
                                  const char *func_name);

void init_uniform_box(SPHSystem *sph, int nx, int ny);

void init_sod_2d(SPHSystem *sph, int nx, int ny);
void init_sph_parameter(SPHSystem *sph);
void calculate_position(int N, double* pos_x, double* pos_y, double x_min, double x_max, double y_min, double y_max, char tail);
void init_sod_2d_2(SPHSystem *sph, double x, double y, double mass);
void init_sod_2d_3(SPHSystem *sph, double x, double y, double mass);
void init_KH(SPHSystem *sph, int nx, int ny);

#endif
