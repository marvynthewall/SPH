# include "sph_system.h"


void compute_forces_dummy(SPHSystem2D *sph)
{
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].ax = 2.0;
        sph->particles[i].ay = 0.0;
    }
}



int main(int argc, char**argv)
{
    printf("========== SPH Integrator Test ==========\n");

    SPHSystem2D sph;
    allocate_sph_system(&sph, 1);


    sph.particles[0].h  = 1.0;
    sph.particles[0].cs = 1.0;

    sph.particles[0].x  = 0.0;
    sph.particles[0].y  = 0.0;

    sph.particles[0].vx = 1.0;
    sph.particles[0].vy = 0.0;

    sph.particles[0].ax = 2.0;
    sph.particles[0].ay = 0.0;


    sph.time  =  0.0;
    sph.t_end = 1.0;
    sph.cfl = 0.1;
    
    if (strcmp(argv[1], "euler") == 0) { printf("Running Euler Method\n"); }
    else if (strcmp(argv[1], "kdk") == 0) { printf("Running KDK Method\n"); }
    
    int step = 0;
    while (sph.time < sph.t_end) {
        printf("step = %2d, time = %.2f, dt = %.2f ",
            step, sph.time, sph.dt);
        printf("x = %.2f, " "vx = %.2f, ax = %.2f\n",
            sph.particles[0].x,
            sph.particles[0].vx,
            sph.particles[0].ax
        );

        if(strcmp(argv[1], "euler") == 0){
            step_euler(&sph);
        }else if(strcmp(argv[1], "kdk") == 0){
            step_leapfrog_kdk(&sph, compute_forces_dummy);
        }
        step++;
    }

    free_sph_system(&sph);

    return 0;
}