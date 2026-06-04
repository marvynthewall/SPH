#include "sph_system.h"

void write_csv(SPHSystem *sph, const char *filename) {
  FILE *fp = fopen(filename, "w");

  if (fp == NULL) {
    fprintf(stderr, "Error: cannot open file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "id,x,y,z,vx,vy,vz,ax,ay,az,m,rho,P,u,h,cs\n");

  for (int i = 0; i < sph->N; i++) {
    fprintf(fp,
            "%d,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%."
            "4e,%.4e,%.4e\n",
            sph->particles[i].id, sph->particles[i].x, sph->particles[i].y,
            sph->particles[i].z, sph->particles[i].vx, sph->particles[i].vy,
            sph->particles[i].vz, sph->particles[i].ax, sph->particles[i].ay,
            sph->particles[i].az, sph->particles[i].mass, sph->particles[i].rho,
            sph->particles[i].pressure, sph->particles[i].u,
            sph->particles[i].h, sph->particles[i].cs);
  }

  fclose(fp);
}

void write_binary(SPHSystem *sph, const char *filename) {
    FILE *fp = fopen(filename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int N = sph->N;
    int ncols = 16;

    double *data = malloc((size_t)N * ncols * sizeof(double));

    if (data == NULL) {
        fprintf(stderr, "Error: memory allocation failed in write_binary.\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        Particle *p = &sph->particles[i];

        data[i * ncols + 0]  = (double)p->id;
        data[i * ncols + 1]  = p->x;
        data[i * ncols + 2]  = p->y;
        data[i * ncols + 3]  = p->z;
        data[i * ncols + 4]  = p->vx;
        data[i * ncols + 5]  = p->vy;
        data[i * ncols + 6]  = p->vz;
        data[i * ncols + 7]  = p->ax;
        data[i * ncols + 8]  = p->ay;
        data[i * ncols + 9]  = p->az;
        data[i * ncols + 10] = p->mass;
        data[i * ncols + 11] = p->rho;
        data[i * ncols + 12] = p->pressure;
        data[i * ncols + 13] = p->u;
        data[i * ncols + 14] = p->h;
        data[i * ncols + 15] = p->cs;
    }

    fwrite(&N, sizeof(int), 1, fp);
    fwrite(&ncols, sizeof(int), 1, fp);
    fwrite(data, sizeof(double), (size_t)N * ncols, fp);

    free(data);
    fclose(fp);
}