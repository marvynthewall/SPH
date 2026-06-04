#include "force.h"
#include "sph_system.h"

// this is for 1-layer loop
// factor, const = 2.0 for 2D, 3.0 for 3D
// Eq(23)
void compute_pressure_soundspeed_factor(SPHSystem *sph) {
  double GG1 = sqrt(sph->gamma * (sph->gamma - 1.0));
  double inv_dim = 1.0 / (double)sph->dim;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].pressure =
        (sph->gamma - 1.0) * sph->particles[i].rho * sph->particles[i].u;
    sph->particles[i].cs = GG1 * sqrt(sph->particles[i].u);
    sph->particles[i].factor =
        1.0 / (1.0 + inv_dim * sph->particles[i].h / sph->particles[i].rho *
                         sph->particles[i].drho_dh);
  }
  return;
}

// make sure that this is always inline
// inline: save time for jumping around different function part.
// static: only inside this file, while compiling, this is better for
// optimization
// __attribute__((always_inline))
// when compiling, compiler can discard this after compiling this part, and
// optimize it directly. inline is only a "suggestion", especially in O3,
// sometimes it changes idea, since this is a long function however in our
// simulation, this need to be inline!
__attribute__((always_inline)) static inline void
compute_pairwise_physics(Particle *p_i, Particle *p_j, SPHSystem *sph) {
  double dx = p_i->x - p_j->x;
  double dy = p_i->y - p_j->y;
  double r = sqrt(dx * dx + dy * dy);
  if (r < 1e-12)
    return; // avoid error

  double max_h = fmax(p_i->h, p_j->h);
  if (r > max_h)
    return; // Viscosity & Force only active within max_h

  double dvx = p_i->vx - p_j->vx;
  double dvy = p_i->vy - p_j->vy;

  // get the W_ij(h_i) and W_ij(h_j)
  double W_i, dWdr_i, dWdh_i;
  cubic_spline_kernel(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
  double W_j, dWdr_j, dWdh_j;
  cubic_spline_kernel(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

  // --- 1. Pressure Force ---
  // Pressure force, Eq(22)
  double term_i = p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
  double term_j = p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;
  double scalar_force = p_j->mass * (term_i + term_j);

  double ax = -scalar_force * (dx / r);
  double ay = -scalar_force * (dy / r);

  p_i->ax += ax;
  p_i->ay += ay;

  // --- 2. Thermal Energy Evolution (Pressure Part) ---
  // update the specific thermal energy evolution, Eq(25)
  double inner_product_v_dW_i =
      dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
  double inner_product_v_dW_j =
      dvx * (dWdr_j * dx / r) +
      dvy * (dWdr_j * dy / r); // 必須保留，avg_inner 會用到

  p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * p_j->mass *
               inner_product_v_dW_i;

  // --- 3. Artificial Viscosity ---
  double r_dot_v = dx * dvx + dy * dvy;
  double Vax = 0.0, Vay = 0.0;

  if (r_dot_v < 0.0) {
    // Eq 31 mu_ij calculation
    double h_ij = (p_i->h + p_j->h) / 2.0;
    double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
    // Eq 30 PI_ij calculation
    double c_ij = (p_i->cs + p_j->cs) / 2.0;
    double rho_ij = (p_i->rho + p_j->rho) / 2.0;
    double PI_ij =
        (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

    // update the specific thermal energy evolution
    // Eq(29)
    // inner product with the avg W is the avg of inner products
    double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
    p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

    // update the acceleration
    // Eq(26)
    Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dx/r);
    Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dy/r);
    p_i->ax += Vax;
    p_i->ay += Vay;
#ifndef _OPENMP
    p_j->dudt += p_i->mass / 2.0 * PI_ij * avg_inner;
#endif
}

  // --- 4. Reaction force only for single core
#ifndef _OPENMP
  double mass_ratio = p_i->mass / p_j->mass;
  p_j->ax -= (ax + Vax) * mass_ratio;
  p_j->ay -= (ay + Vay) * mass_ratio;
  p_j->dudt += p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * p_i->mass *
               inner_product_v_dW_j;
#endif
}

void compute_force(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  // initialize the ax,ay,dudt in every particle
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }
#ifdef _OPENMP
#pragma omp parallel for
  // openmp, run every i to j, for parallelization
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    for (int j = 0; j < sph->N; j++) {
      if (i == j)
        continue;
      Particle *p_j = &sph->particles[j];
      compute_pairwise_physics(p_i, p_j, sph);
    }
  }
#else
  // run only pairs (half the number of the computation)
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    for (int j = i + 1; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j];
      compute_pairwise_physics(p_i, p_j, sph);
    }
  }
#endif
}


// Computing force, symmetric 2-layers loop with boundary conditions
void compute_force_xreflective_yperiodic_celllist(SPHSystem *sph) {
  // 1. Initialize acceleration and thermal energy evolution rate for all real particles
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

  // 2. Double full-loop: For OpenMP thread-safety, each thread only writes to p_i.
  // The workload is distributed and parallelized across all CPU cores.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];

    // current cell
    int cx_i = (int)(p_i->x / sph->cell_size);
    int cy_i = (int)(p_i->y / sph->cell_size);

// --- 遍歷相鄰的 9 個網格 (包含自己所在的網格) ---
    for (int d_cy = -1; d_cy <= 1; d_cy++) {
      for (int d_cx = -1; d_cx <= 1; d_cx++) {
        
        int cx = cx_i + d_cx;
        int cy = cy_i + d_cy;

        // 【Y 方向週期性邊界】網格座標折返
        if (cy < 0) cy += sph->num_cells_y;
        else if (cy >= sph->num_cells_y) cy -= sph->num_cells_y;

        // 【X 方向反射邊界】超出邊界代表是牆壁，沒有真實網格，直接跳過
        if (cx < 0 || cx >= sph->num_cells_x) continue;

        int cell_index = cx + cy * sph->num_cells_x;
        int j = sph->head[cell_index];

        // --- 走訪該網格內的 Linked List ---
        while (j != -1) {
          Particle *p_j = &sph->particles[j];

          // Y-Periodic 距離修正 (Minimum Image Convention)
          double dy = p_i->y - p_j->y;
          if (dy > 0.5 * sph->box_size_y) dy -= sph->box_size_y;
          else if (dy < -0.5 * sph->box_size_y) dy += sph->box_size_y;

          // ==========================================
          // 形態 0：真實粒子交互作用
          // ==========================================
          if (i != j) { // 排除自己對自己產生真實作用力
            Particle ghost_j = *p_j;
            ghost_j.y  = p_i->y - dy; // 修正 Y 座標以符合週期距離
            // x, vx, vy 皆維持真實狀態
            compute_pairwise_physics(p_i, &ghost_j, sph);
          }

          // ==========================================
          // 形態 1：左反射牆鏡像 (Left Mirror)
          // ==========================================
          // 物理防呆剪枝：只有當 i 靠近左邊界時才計算。
          // ⚠️ 這裡必須用全局的 cell_size (即 max_h)，而非 p_i->h，
          // 以免 p_i 的 h 很小，但 p_j 的 h 很大時，錯失了交互作用。
          if (p_i->x < sph->cell_size) {
            Particle ghost_j = *p_j;
            ghost_j.x  = -p_j->x;
            ghost_j.y  = p_i->y - dy;
            ghost_j.vx = -p_j->vx; // 撞牆反彈，X 速度反轉
            // vy 維持不變
            compute_pairwise_physics(p_i, &ghost_j, sph);
          }

          // ==========================================
          // 形態 2：右反射牆鏡像 (Right Mirror)
          // ==========================================
          if (p_i->x > sph->box_size_x - sph->cell_size) {
            Particle ghost_j = *p_j;
            ghost_j.x  = 2.0 * sph->box_size_x - p_j->x;
            ghost_j.y  = p_i->y - dy;
            ghost_j.vx = -p_j->vx; // 撞牆反彈，X 速度反轉
            compute_pairwise_physics(p_i, &ghost_j, sph);
          }

          j = sph->next[j]; 
        } // end while (j != -1)
      }
    } // end for d_cx, d_cy
  }

  // 3. Handle Self-mirror Interactions
  // When a particle gets close to the left or right walls, it interacts with its own mirror image.
  // Each particle calculates this independently, so it can be parallelized with OpenMP.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    
    // Two self-mirror positions: left-wall mirror and right-wall mirror
    double dxs[2] = {2.0 * p_i->x, 2.0 * (p_i->x - sph->box_size_x)};
    double dvxs[2] = {2.0 * p_i->vx, 2.0 * p_i->vx};

    for (int k = 0; k < 2; k++) {
      // Create a ghost particle representing the particle's own mirror image
      Particle ghost_self = *p_i;
      
      ghost_self.x  = p_i->x - dxs[k];
      ghost_self.vx = p_i->vx - dvxs[k];
      // Y-direction perfectly overlaps with itself; relative distance and velocity are 0
      ghost_self.y  = p_i->y;
      ghost_self.vy = p_i->vy; 

      // Invoke the same core physics function
      compute_pairwise_physics(p_i, &ghost_self, sph);
    }
  }
}


// Computing force, symmetric 2-layers loop with boundary conditions
void compute_force_xreflective_yperiodic(SPHSystem *sph) {
  // 1. Initialize acceleration and thermal energy evolution rate for all real particles
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

  // 2. Double full-loop: For OpenMP thread-safety, each thread only writes to p_i.
  // The workload is distributed and parallelized across all CPU cores.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {
      // Skip self-interaction; self-mirror interaction is handled separately below
      if (i == j) continue; 
      
      Particle *p_j = &sph->particles[j];

      // --- Handle Y-direction Periodic Boundary Condition ---
      // Calculate relative y-displacement using the Minimum Image Convention
      double dy = p_i->y - p_j->y;
      dy -= (dy >  0.5 * sph->box_size_y) ? sph->box_size_y : 0.0;
      dy += (dy < -0.5 * sph->box_size_y) ? sph->box_size_y : 0.0;

      // --- Handle X-direction Reflective/Mirror Boundary Condition ---
      // Define the 3 potential spatial configurations of particle j relative to i:
      // Index 0: Original relative position
      // Index 1: Mirror image behind the left boundary wall
      // Index 2: Mirror image behind the right boundary wall
      double dxs[3] = {
          p_i->x - p_j->x,                         
          p_i->x + p_j->x,                         
          p_i->x + p_j->x - 2.0 * sph->box_size_x  
      };
      
      // Mirror particles must have their X-velocity inverted to correctly simulate 
      // the relative rebound velocity against the wall: vx_i - (-vx_j) = vx_i + vx_j
      double dvxs[3] = {
          p_i->vx - p_j->vx, 
          p_i->vx + p_j->vx, 
          p_i->vx + p_j->vx  
      };

      // Traverse all 3 spatial configurations
      for (int k = 0; k < 3; k++) {
        // Create a temporary ghost particle j by copying the properties of real particle j
        Particle ghost_j = *p_j;

        // Reconstruct the absolute coordinates and velocities of ghost_j based on geometric relations
        ghost_j.x  = p_i->x - dxs[k];
        ghost_j.y  = p_i->y - dy;
        ghost_j.vx = p_i->vx - dvxs[k];
        ghost_j.vy = p_i->vy - (p_i->vy - p_j->vy); // Maintain original dvy

        // Reuse the core pairwise physics function.
        // Note: When _OPENMP is defined, this function ONLY updates the first argument (p_i).
        // The code wrapped inside #ifndef _OPENMP will be automatically disabled, 
        // preventing any race conditions or unnecessary writes to ghost_j.
        compute_pairwise_physics(p_i, &ghost_j, sph);
      }
    }
  }

  // 3. Handle Self-mirror Interactions
  // When a particle gets close to the left or right walls, it interacts with its own mirror image.
  // Each particle calculates this independently, so it can be parallelized with OpenMP.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    
    // Two self-mirror positions: left-wall mirror and right-wall mirror
    double dxs[2] = {2.0 * p_i->x, 2.0 * (p_i->x - sph->box_size_x)};
    double dvxs[2] = {2.0 * p_i->vx, 2.0 * p_i->vx};

    for (int k = 0; k < 2; k++) {
      // Create a ghost particle representing the particle's own mirror image
      Particle ghost_self = *p_i;
      
      ghost_self.x  = p_i->x - dxs[k];
      ghost_self.vx = p_i->vx - dvxs[k];
      // Y-direction perfectly overlaps with itself; relative distance and velocity are 0
      ghost_self.y  = p_i->y;
      ghost_self.vy = p_i->vy; 

      // Invoke the same core physics function
      compute_pairwise_physics(p_i, &ghost_self, sph);
    }
  }
}


void compute_force_xperiodic_yperiodic(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sph->N; i++) {
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].dudt = 0.0;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];

    for (int j = 0; j < sph->N; j++) {
      if (i == j)
        continue;
      Particle *p_j = &sph->particles[j];

      double dx = p_i->x - p_j->x;
      if (dx > 0.5 * sph->box_size_x)
        dx -= sph->box_size_x;
      else if (dx < -0.5 * sph->box_size_x)
        dx += sph->box_size_x;

      double dy = p_i->y - p_j->y;
      if (dy > 0.5 * sph->box_size_y)
        dy -= sph->box_size_y;
      else if (dy < -0.5 * sph->box_size_y)
        dy += sph->box_size_y;

      double r = sqrt(dx * dx + dy * dy);
      double max_h = fmax(p_i->h, p_j->h);

      if (r < 1e-12 || r > max_h)
        continue;

      double dvx = p_i->vx - p_j->vx;
      double dvy = p_i->vy - p_j->vy;

      double W_i, dWdr_i, dWdh_i;
      cubic_spline_kernel(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
      double W_j, dWdr_j, dWdh_j;
      cubic_spline_kernel(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

      double term_i =
          p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
      double term_j =
          p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;

      double scalar_force = p_j->mass * (term_i + term_j);

      double ax = -scalar_force * (dx / r);
      double ay = -scalar_force * (dy / r);

      p_i->ax += ax;
      p_i->ay += ay;

      double inner_product_v_dW_i =
          dvx * (dWdr_i * dx / r) + dvy * (dWdr_i * dy / r);
      double inner_product_v_dW_j =
          dvx * (dWdr_j * dx / r) + dvy * (dWdr_j * dy / r);
      p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) *
                   p_j->mass * inner_product_v_dW_i;

      double r_dot_v = dx * dvx + dy * dvy;
      if (r_dot_v < 0.0) {
        double h_ij = (p_i->h + p_j->h) / 2.0;
        double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
        double c_ij = (p_i->cs + p_j->cs) / 2.0;
        double rho_ij = (p_i->rho + p_j->rho) / 2.0;
        double PI_ij =
            (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

        double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
        p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

        double Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dx / r);
        double Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j) / 2.0 * dy / r);
        p_i->ax += Vax;
        p_i->ay += Vay;
      }
    }
  }
}

__attribute__((always_inline)) static inline void compute_pairwise_physics_3d(
      Particle * p_i, Particle * p_j, SPHSystem * sph) {
    double dx = p_i->x - p_j->x;
    double dy = p_i->y - p_j->y;
    double dz = p_i->z - p_j->z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-12)
      return; // avoid error

    double max_h = fmax(p_i->h, p_j->h);
    if (r > max_h)
      return; // Viscosity & Force only active within max_h

    double dvx = p_i->vx - p_j->vx;
    double dvy = p_i->vy - p_j->vy;
    double dvz = p_i->vz - p_j->vz;

    // get the W_ij(h_i) and W_ij(h_j)
    double W_i, dWdr_i, dWdh_i;
    cubic_spline_kernel_3d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);
    double W_j, dWdr_j, dWdh_j;
    cubic_spline_kernel_3d(r, p_j->h, &W_j, &dWdr_j, &dWdh_j);

    // --- 1. Pressure Force ---
    // Pressure force, Eq(22)
    double term_i =
        p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) * dWdr_i;
    double term_j =
        p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) * dWdr_j;
    double scalar_force = p_j->mass * (term_i + term_j);

    double ax = -scalar_force * (dx / r);
    double ay = -scalar_force * (dy / r);
    double az = -scalar_force * (dz / r);

    p_i->ax += ax;
    p_i->ay += ay;
    p_i->az += az;

    // --- 2. Thermal Energy Evolution (Pressure Part) ---
    // update the specific thermal energy evolution, Eq(25)
    double inner_product_v_dW_i = dvx * (dWdr_i * dx / r) +
                                  dvy * (dWdr_i * dy / r) +
                                  dvz * (dWdr_i * dz / r);
    double inner_product_v_dW_j = dvx * (dWdr_j * dx / r) +
                                  dvy * (dWdr_j * dy / r) +
                                  dvz * (dWdr_j * dz / r);

    p_i->dudt += p_i->factor * p_i->pressure / (p_i->rho * p_i->rho) *
                 p_j->mass * inner_product_v_dW_i;

    // --- 3. Artificial Viscosity ---
    double r_dot_v = dx * dvx + dy * dvy + dz * dvz;
    double Vax = 0.0, Vay = 0.0, Vaz = 0.0;

    if (r_dot_v < 0.0) {
      // Eq 31 mu_ij calculation
      double h_ij = (p_i->h + p_j->h) / 2.0;
      double mu_ij = h_ij * r_dot_v / (r * r + sph->epsilon * (h_ij * h_ij));
      // Eq 30 PI_ij calculation
      double c_ij = (p_i->cs + p_j->cs) / 2.0;
      double rho_ij = (p_i->rho + p_j->rho) / 2.0;
      double PI_ij =
          (-sph->alpha * c_ij * mu_ij + sph->beta * mu_ij * mu_ij) / rho_ij;

      // update the specific thermal energy evolution
      // Eq(29)
      // inner product with the avg W is the avg of inner products
      double avg_inner = (inner_product_v_dW_i + inner_product_v_dW_j) / 2.0;
      p_i->dudt += p_j->mass / 2.0 * PI_ij * avg_inner;

      // update the acceleration
      // Eq(26)
      Vax = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dx/r);
      Vay = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dy/r);
      Vaz = -p_j->mass * PI_ij * ((dWdr_i + dWdr_j)/2.0 * dz/r);
      p_i->ax += Vax;
      p_i->ay += Vay;
      p_i->az += Vaz;
#ifndef _OPENMP
      p_j->dudt += p_i->mass / 2.0 * PI_ij * avg_inner;
#endif
    }

    // --- 4. Reaction force only for single core
#ifndef _OPENMP
    double mass_ratio = p_i->mass / p_j->mass;
    p_j->ax -= (ax + Vax) * mass_ratio;
    p_j->ay -= (ay + Vay) * mass_ratio;
    p_j->az -= (az + Vaz) * mass_ratio;
    p_j->dudt += p_j->factor * p_j->pressure / (p_j->rho * p_j->rho) *
                 p_i->mass * inner_product_v_dW_j;
#endif
}

void compute_force_3d(SPHSystem *sph) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  // initialize the ax,ay,dudt in every particle
  for (int i = 0; i < sph->N; i++){
    sph->particles[i].ax = 0.0;
    sph->particles[i].ay = 0.0;
    sph->particles[i].az = 0.0;
    sph->particles[i].dudt = 0.0;
  }
#ifdef _OPENMP
#pragma omp parallel for
  // openmp, run every i to j, for parallelization
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    for (int j = 0; j < sph->N; j++) {
      if (i == j)
        continue;
      Particle *p_j = &sph->particles[j];
      compute_pairwise_physics_3d(p_i, p_j, sph);
    }
  }
#else
  // run only pairs (half the number of the computation)
  for (int i = 0; i < sph->N; i++) {
    Particle *p_i = &sph->particles[i];
    for (int j = i + 1; j < sph->N; j++) {
      Particle *p_j = &sph->particles[j];
      compute_pairwise_physics_3d(p_i, p_j, sph);
    }
  }
#endif
}


void compute_force_xreflective_yperiodic_zperiodic_3d(SPHSystem *sph)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].ax = 0.0;
        sph->particles[i].ay = 0.0;
        sph->particles[i].az = 0.0;
        sph->particles[i].dudt = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        for (int j = 0; j < sph->N; j++) {
            if (i == j) continue;
            

            Particle *p_j = &sph->particles[j];

            double dx = p_i->x - p_j->x;

            double dy = p_i->y - p_j->y;
            if (dy > 0.5 * sph->box_size_y) {
                dy -= sph->box_size_y;
            } else if (dy < -0.5 * sph->box_size_y) {
                dy += sph->box_size_y;
            }

            double dz = p_i->z - p_j->z;
            if (dz > 0.5 * sph->box_size_z) {
                dz -= sph->box_size_z;
            } else if (dz < -0.5 * sph->box_size_z) {
                dz += sph->box_size_z;
            }

            double r = sqrt(dx * dx + dy * dy + dz * dz);

            double max_h = fmax(p_i->h, p_j->h);

            if (r < 1.0e-12 || r > max_h) {
                continue;
            }

            double W_i = 0.0;
            double dWdr_i = 0.0;
            double dWdh_i = 0.0;

            cubic_spline_kernel_3d(r, p_i->h, &W_i, &dWdr_i, &dWdh_i);

            double gradWix = dWdr_i * dx / r;
            double gradWiy = dWdr_i * dy / r;
            double gradWiz = dWdr_i * dz / r;

            double dvx = p_i->vx - p_j->vx;
            double dvy = p_i->vy - p_j->vy;
            double dvz = p_i->vz - p_j->vz;

            double pressure_term =
                p_i->factor * p_i->pressure / (p_i->rho * p_i->rho)
              + p_j->factor * p_j->pressure / (p_j->rho * p_j->rho);

            p_i->ax += -p_j->mass * pressure_term * gradWix;
            p_i->ay += -p_j->mass * pressure_term * gradWiy;
            p_i->az += -p_j->mass * pressure_term * gradWiz;

            double inner_product_v_dW_i =
                dvx * gradWix
              + dvy * gradWiy
              + dvz * gradWiz;

            p_i->dudt +=
                p_i->factor
              * p_i->pressure
              / (p_i->rho * p_i->rho)
              * p_j->mass
              * inner_product_v_dW_i;

            double r_dot_v =
                dx * dvx
              + dy * dvy
              + dz * dvz;

            if (r_dot_v < 0.0) {

                double h_ij = 0.5 * (p_i->h + p_j->h);

                double mu_ij =
                    h_ij * r_dot_v
                  / (r * r + sph->epsilon * h_ij * h_ij);

                double c_ij = 0.5 * (p_i->cs + p_j->cs);
                double rho_ij = 0.5 * (p_i->rho + p_j->rho);

                double PI_ij =
                    (-sph->alpha * c_ij * mu_ij
                     + sph->beta * mu_ij * mu_ij)
                  / rho_ij;

                p_i->dudt +=
                    0.5
                  * p_j->mass
                  * PI_ij
                  * inner_product_v_dW_i;

                p_i->ax += -p_j->mass * PI_ij * gradWix;
                p_i->ay += -p_j->mass * PI_ij * gradWiy;
                p_i->az += -p_j->mass * PI_ij * gradWiz;
            }
        }
    }
}




void compute_force_xreflective_yzperiodic_celllist_3d(SPHSystem *sph)
{
    /*
     * Precondition:
     * build_cell_list_3d(sph) should already be called before this function.
     *
     * In your integrator, this is usually done inside:
     * compute_density_xreflective_yzperiodic_celllist_3d(sph)
     *
     * If you want this force function to be fully independent,
     * you can uncomment the following line:
     *
     * build_cell_list_3d(sph);
     */

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < sph->N; i++) {
        sph->particles[i].ax = 0.0;
        sph->particles[i].ay = 0.0;
        sph->particles[i].az = 0.0;
        sph->particles[i].dudt = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        int cx_i = (int)(p_i->x / sph->cell_size);
        int cy_i = (int)(p_i->y / sph->cell_size);
        int cz_i = (int)(p_i->z / sph->cell_size);

        if (cx_i < 0) cx_i = 0;
        if (cx_i >= sph->num_cells_x) cx_i = sph->num_cells_x - 1;

        if (cy_i < 0) cy_i = 0;
        if (cy_i >= sph->num_cells_y) cy_i = sph->num_cells_y - 1;

        if (cz_i < 0) cz_i = 0;
        if (cz_i >= sph->num_cells_z) cz_i = sph->num_cells_z - 1;

        /*
         * Search 27 neighboring cells:
         * dx direction: reflective wall, skip outside cells
         * y,z directions: periodic wrapping
         */
        for (int d_cz = -1; d_cz <= 1; d_cz++) {
            for (int d_cy = -1; d_cy <= 1; d_cy++) {
                for (int d_cx = -1; d_cx <= 1; d_cx++) {

                    int cx = cx_i + d_cx;
                    int cy = cy_i + d_cy;
                    int cz = cz_i + d_cz;

                    // Y-periodic cell wrapping
                    if (cy < 0) {
                        cy += sph->num_cells_y;
                    } else if (cy >= sph->num_cells_y) {
                        cy -= sph->num_cells_y;
                    }

                    // Z-periodic cell wrapping
                    if (cz < 0) {
                        cz += sph->num_cells_z;
                    } else if (cz >= sph->num_cells_z) {
                        cz -= sph->num_cells_z;
                    }

                    // X-reflective: no real cells outside x walls
                    if (cx < 0 || cx >= sph->num_cells_x) {
                        continue;
                    }

                    int cell_index =
                        cx
                      + cy * sph->num_cells_x
                      + cz * sph->num_cells_x * sph->num_cells_y;

                    int j = sph->head[cell_index];

                    while (j != -1) {

                        Particle *p_j = &sph->particles[j];

                        /*
                         * Y/Z periodic minimum image.
                         * Then reconstruct ghost_j.y and ghost_j.z
                         * so that compute_pairwise_physics_3d() sees
                         * the corrected displacement.
                         */
                        double dy = p_i->y - p_j->y;
                        if (dy > 0.5 * sph->box_size_y) {
                            dy -= sph->box_size_y;
                        } else if (dy < -0.5 * sph->box_size_y) {
                            dy += sph->box_size_y;
                        }

                        double dz = p_i->z - p_j->z;
                        if (dz > 0.5 * sph->box_size_z) {
                            dz -= sph->box_size_z;
                        } else if (dz < -0.5 * sph->box_size_z) {
                            dz += sph->box_size_z;
                        }

                        // ==================================================
                        // 1. Real particle interaction
                        // ==================================================
                        if (i != j) {
                            Particle ghost_j = *p_j;

                            ghost_j.x = p_j->x;
                            ghost_j.y = p_i->y - dy;
                            ghost_j.z = p_i->z - dz;

                            /*
                             * Velocities are unchanged for periodic y/z.
                             */
                            ghost_j.vx = p_j->vx;
                            ghost_j.vy = p_j->vy;
                            ghost_j.vz = p_j->vz;

                            compute_pairwise_physics_3d(p_i, &ghost_j, sph);
                        }

                        // ==================================================
                        // 2. Left reflective mirror across x = 0
                        // x_j,ghost = -x_j
                        // vx_j,ghost = -vx_j
                        // ==================================================
                        if (p_i->x < sph->cell_size) {
                            Particle ghost_j = *p_j;

                            ghost_j.x = -p_j->x;
                            ghost_j.y = p_i->y - dy;
                            ghost_j.z = p_i->z - dz;

                            ghost_j.vx = -p_j->vx;
                            ghost_j.vy =  p_j->vy;
                            ghost_j.vz =  p_j->vz;

                            compute_pairwise_physics_3d(p_i, &ghost_j, sph);
                        }

                        // ==================================================
                        // 3. Right reflective mirror across x = box_size_x
                        // x_j,ghost = 2Lx - x_j
                        // vx_j,ghost = -vx_j
                        // ==================================================
                        if (p_i->x > sph->box_size_x - sph->cell_size) {
                            Particle ghost_j = *p_j;

                            ghost_j.x = 2.0 * sph->box_size_x - p_j->x;
                            ghost_j.y = p_i->y - dy;
                            ghost_j.z = p_i->z - dz;

                            ghost_j.vx = -p_j->vx;
                            ghost_j.vy =  p_j->vy;
                            ghost_j.vz =  p_j->vz;

                            compute_pairwise_physics_3d(p_i, &ghost_j, sph);
                        }

                        j = sph->next[j];
                    }
                }
            }
        }
    }

    /*
     * 4. Self-mirror interactions.
     *
     * This handles particle i interacting with its own reflected image
     * near the left and right x-walls.
     */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < sph->N; i++) {

        Particle *p_i = &sph->particles[i];

        // Left self-mirror
        if (p_i->x < sph->cell_size) {
            Particle ghost_self = *p_i;

            ghost_self.x  = -p_i->x;
            ghost_self.y  =  p_i->y;
            ghost_self.z  =  p_i->z;

            ghost_self.vx = -p_i->vx;
            ghost_self.vy =  p_i->vy;
            ghost_self.vz =  p_i->vz;

            compute_pairwise_physics_3d(p_i, &ghost_self, sph);
        }

        // Right self-mirror
        if (p_i->x > sph->box_size_x - sph->cell_size) {
            Particle ghost_self = *p_i;

            ghost_self.x  = 2.0 * sph->box_size_x - p_i->x;
            ghost_self.y  = p_i->y;
            ghost_self.z  = p_i->z;

            ghost_self.vx = -p_i->vx;
            ghost_self.vy =  p_i->vy;
            ghost_self.vz =  p_i->vz;

            compute_pairwise_physics_3d(p_i, &ghost_self, sph);
        }
    }
}