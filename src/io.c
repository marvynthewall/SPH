#include "sph_system.h"

void write_csv(SPHSystem *sph, const char *filename)
{
    FILE *fp = fopen(filename, "w");

    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "id,x,y,z,vx,vy,vz,ax,ay,az,m,rho,P,u,h,cs\n");

    for (int i = 0; i < sph->N; i++) {
        fprintf(fp,
                "%d,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                sph->particles[i].id,
                sph->particles[i].x,
                sph->particles[i].y,
                sph->particles[i].z,
                sph->particles[i].vx,
                sph->particles[i].vy,
                sph->particles[i].vz,
                sph->particles[i].ax,
                sph->particles[i].ay,
                sph->particles[i].az,
                sph->particles[i].mass,
                sph->particles[i].rho,
                sph->particles[i].pressure,
                sph->particles[i].u,
                sph->particles[i].h,
                sph->particles[i].cs);
    }

    fclose(fp);
}
