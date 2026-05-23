#include "sph_system.h"

int main() {
  printf("========== SPH Density Test ==========\n");

  Particle particles[2];

  // 初始化粒子 0 (放在原點)
  particles[0].id = 0;
  particles[0].x = 0.0;
  particles[0].y = 0.0;
  particles[0].mass = 1.0;
  particles[0].h = 1.0;
  particles[0].density = 0.0;

  // 初始化粒子 1 (放在 x=0.5，跟粒子0距離 0.5)
  particles[1].id = 1;
  particles[1].x = 0.5;
  particles[1].y = 0.0;
  particles[1].mass = 1.0;
  particles[1].h = 1.0;
  particles[1].density = 0.0;

  // 呼叫密度計算函數
  compute_density(particles, 2);

  // 印出計算結果
  printf("Particle 0 Density: %f\n", particles[0].density);
  printf("Particle 1 Density: %f\n", particles[1].density);
  printf("Note: 它們的密度應該要一樣，因為質量和相對距離相同！\n");

  printf("\n================ Test Finished ================\n");
  return 0;
}
