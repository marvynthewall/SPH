# include "sph_system.h"

void write_csv(SPHSystem2D * sph, const char * filename){
    FILE * fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "id,x,y,vx,vy,ax,ay,m,rho,P,u,h\n");

    for (int i = 0; i < sph->N; i++) {
        fprintf(fp,
                "%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                i,
                sph->x[i],
                sph->y[i],
                sph->vx[i],
                sph->vy[i],
                sph->ax[i],
                sph->ay[i],
                sph->m[i],
                sph->rho[i],
                sph->P[i],
                sph->u[i],
                sph->h[i]);
    }

    fclose(fp);
}