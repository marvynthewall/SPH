#include "sph_system.h"
#include <math.h>


void test_pressure_soundspeed() {
    SPHSystem sph;
    allocate_sph_system(&sph, 1);
    sph.particles[0].rho = 1.0;
    sph.particles[0].u = 1.0;
    sph.particles[0].h = 0.1;
    sph.particles[0].drho_dh = 0.0;

    compute_pressure_soundspeed_factor(&sph);

    // γ=2: P = (2-1)*1.0*1.0 = 1.0
    assert(fabs(sph.particles[0].pressure - 1.0) < 1e-10);
    
    // cs = sqrt(γ*(γ-1)*u) = sqrt(2.0*1.0*1.0) = sqrt(2)
    double expected_cs = sqrt(GAMMA * (GAMMA - 1) * 1.0);
    assert(fabs(sph.particles[0].cs - expected_cs) < 1e-10);

    // factor = 1.0
    assert(fabs(sph.particles[0].factor - 1.0) < 1e-10);
    
    printf("PASS: pressure_soundspeed\n");
}

void test_force_two_particle_head_on_collision(void){
    SPHSystem sph;
    allocate_sph_system(&sph, 2);
    
    for (int i = 0; i < 2; i++) {
        sph.particles[i].mass = 1.0;
        sph.particles[i].rho = 1.0;
        sph.particles[i].pressure = 1.0;
        sph.particles[i].u = 1.0;
        sph.particles[i].h = 1.0;
        sph.particles[i].cs = 1.0;
        sph.particles[i].factor = 1.0;
        sph.particles[i].drho_dh = 0.0;
    }

    sph.particles[0].x = 0.0;  sph.particles[0].y = 0.0;
    sph.particles[1].x = 0.5;  sph.particles[1].y = 0.0;

    // 迎面相撞（相對速度 dvx = 1.0 - (-1.0) = 2.0）
    sph.particles[0].vx = 1.0;   sph.particles[0].vy = 0.0;
    sph.particles[1].vx = -1.0;  sph.particles[1].vy = 0.0;

    compute_force(&sph);

    // 6. 印出結果進行驗證
    printf("========== SPH FORCE VERIFICATION ==========\n");
    printf("Particle 0 ax   : %f (預期應為負值，因為被排斥且減速)\n", sph.particles[0].ax);
    printf("Particle 1 ax   : %f (預期應與 P0 大小相同、符號相反)\n", sph.particles[1].ax);
    printf("Particle 0 ay   : %f (預期應精確為 0.0)\n", sph.particles[0].ay);
    printf("Particle 1 ay   : %f (預期應精確為 0.0)\n", sph.particles[1].ay);
    printf("Particle 0 dudt : %f (預期應為正值，因為摩擦生熱)\n", sph.particles[0].dudt);
    printf("Particle 1 dudt : %f (預期應與 P0 完全相同)\n", sph.particles[1].dudt);

    free_sph_system(&sph);
}

int main(void) {
    test_pressure_soundspeed();
    test_force_two_particle_head_on_collision();
   return 0;
}
