#include "sph_system.h"



void test_timestep_single_particle(void)
{
    SPHSystem2D sph;
    allocate_sph_system(&sph, 1);

    sph.particles[0].h  = 1.0;
    sph.particles[0].cs = 1.0;
    sph.particles[0].x  = 0.0;
    sph.particles[0].y  = 0.0;
    sph.particles[0].vx = 1.0;
    sph.particles[0].vy = 0.0;

    sph.cfl = 0.1;

    double dt_basic  = compute_timestep(&sph);
    double dt_signal = compute_timestep_signal_velocity(&sph);

    printf("[N=1]\n");
    printf("dt_basic  = %.15e, expected = %.15e\n", dt_basic,  0.1);
    printf("dt_signal = %.15e, expected = %.15e\n", dt_signal, 0.1);

    free_sph_system(&sph);
}

void test_timestep_two_particles_static(void)
{
    SPHSystem2D sph;
    allocate_sph_system(&sph, 2);

    sph.cfl = 0.1;

    sph.particles[0].x = 0.0;
    sph.particles[0].y = 0.0;
    sph.particles[0].vx = 0.0;
    sph.particles[0].vy = 0.0;
    sph.particles[0].h = 1.0;
    sph.particles[0].cs = 1.0;

    sph.particles[1].x = 1.0;
    sph.particles[1].y = 0.0;
    sph.particles[1].vx = 0.0;
    sph.particles[1].vy = 0.0;
    sph.particles[1].h = 1.0;
    sph.particles[1].cs = 1.0;

    double dt_signal = compute_timestep_signal_velocity(&sph);

    printf("[N=2 static]\n");
    printf("dt_signal = %.15e, expected = %.15e\n", dt_signal, 0.05);

    free_sph_system(&sph);
}

void test_timestep_two_particles_approaching(void)
{
    SPHSystem2D sph;
    allocate_sph_system(&sph, 2);

    sph.cfl = 0.1;

    sph.particles[0].x = 0.0;
    sph.particles[0].y = 0.0;
    sph.particles[0].vx = 1.0;
    sph.particles[0].vy = 0.0;
    sph.particles[0].h = 1.0;
    sph.particles[0].cs = 1.0;

    sph.particles[1].x = 1.0;
    sph.particles[1].y = 0.0;
    sph.particles[1].vx = -1.0;
    sph.particles[1].vy = 0.0;
    sph.particles[1].h = 1.0;
    sph.particles[1].cs = 1.0;

    double dt_signal = compute_timestep_signal_velocity(&sph);

    printf("[N=2 approaching]\n");
    printf("dt_signal = %.15e, expected = %.15e\n", dt_signal, 0.0125);

    free_sph_system(&sph);
}



void init_single_particle_test(SPHSystem2D *sph)
{
    allocate_sph_system(sph, 1);

    sph->particles[0].h  = 1.0;
    sph->particles[0].cs = 1.0;

    sph->particles[0].x  = 0.0;
    sph->particles[0].y  = 0.0;

    sph->particles[0].vx = 1.0;
    sph->particles[0].vy = 0.0;

    sph->particles[0].ax = 0.0;
    sph->particles[0].ay = 0.0;
    sph->particles[0].mass = 1.0;
    sph->particles[0].density = 1.0;
    sph->particles[0].pressure = 0.0;
    sph->particles[0].u = 1.0;

    sph->time  = 0.0;
    sph->dt    = 0.0;
    sph->t_end = 1.0;
    sph->cfl   = 0.1;
}


void compute_forces_dummy(SPHSystem2D *sph)
{
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].ax = 2.0;
        sph->particles[i].ay = 0.0;
    }
}




int main(int argc, char **argv)
{
    printf("========== 1 SPH Timestep Tests ==========\n");

    test_timestep_single_particle();
    test_timestep_two_particles_static();
    test_timestep_two_particles_approaching();

    printf("========== 2 SPH Integrator Unit Tests ==========\n");

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
            step_euler(&sph, compute_timestep);
        }else if(strcmp(argv[1], "kdk") == 0){
            step_leapfrog_kdk(&sph, compute_timestep, compute_forces_dummy);
        }
        step++;
    }

    printf("\nAll integrator tests passed.\n");

    return 0;
}