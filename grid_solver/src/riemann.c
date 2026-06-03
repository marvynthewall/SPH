#include "../include/riemann.h"

void hllc_flux(PrimVar *W_L, PrimVar *W_R, ConsVar *F, double gamma, int dir) {
    double u_L, v_L, u_R, v_R;
    if (dir == 0) { // x-direction
        u_L = W_L->vx; v_L = W_L->vy;
        u_R = W_R->vx; v_R = W_R->vy;
    } else { // y-direction
        u_L = W_L->vy; v_L = W_L->vx;
        u_R = W_R->vy; v_R = W_R->vx;
    }

    double rho_L = W_L->rho;
    double P_L = W_L->P;
    double rho_R = W_R->rho;
    double P_R = W_R->P;

    double a_L = get_soundspeed(gamma, P_L, rho_L);
    double a_R = get_soundspeed(gamma, P_R, rho_R);

    // Compute Roe averages
    double rrho_L = sqrt(rho_L);
    double rrho_R = sqrt(rho_R);
    double u_tilde = (rrho_L * u_L + rrho_R * u_R) / (rrho_L + rrho_R);
    double v_tilde = (rrho_L * v_L + rrho_R * v_R) / (rrho_L + rrho_R);
    
    double H_L = (P_L / (gamma - 1.0) + 0.5 * rho_L * (u_L*u_L + v_L*v_L) + P_L) / rho_L;
    double H_R = (P_R / (gamma - 1.0) + 0.5 * rho_R * (u_R*u_R + v_R*v_R) + P_R) / rho_R;
    double H_tilde = (rrho_L * H_L + rrho_R * H_R) / (rrho_L + rrho_R);
    double a_tilde = sqrt((gamma - 1.0) * (H_tilde - 0.5 * (u_tilde*u_tilde + v_tilde*v_tilde)));

    // Wave speed estimates
    double S_L = fmin(u_L - a_L, u_tilde - a_tilde);
    double S_R = fmax(u_R + a_R, u_tilde + a_tilde);

    double E_L = P_L / (gamma - 1.0) + 0.5 * rho_L * (u_L*u_L + v_L*v_L);
    double F_L_rho = rho_L * u_L;
    double F_L_u = rho_L * u_L * u_L + P_L;
    double F_L_v = rho_L * u_L * v_L;
    double F_L_E = u_L * (E_L + P_L);

    double E_R = P_R / (gamma - 1.0) + 0.5 * rho_R * (u_R*u_R + v_R*v_R);
    double F_R_rho = rho_R * u_R;
    double F_R_u = rho_R * u_R * u_R + P_R;
    double F_R_v = rho_R * u_R * v_R;
    double F_R_E = u_R * (E_R + P_R);

    if (S_L >= 0.0) {
        F->rho = F_L_rho;
        if (dir == 0) { F->px = F_L_u; F->py = F_L_v; }
        else { F->py = F_L_u; F->px = F_L_v; }
        F->E = F_L_E;
        return;
    }
    if (S_R <= 0.0) {
        F->rho = F_R_rho;
        if (dir == 0) { F->px = F_R_u; F->py = F_R_v; }
        else { F->py = F_R_u; F->px = F_R_v; }
        F->E = F_R_E;
        return;
    }

    double S_star = (P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) / 
                    (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

    double F_star_rho, F_star_u, F_star_v, F_star_E;
    if (S_star >= 0.0) {
        double coef = rho_L * (S_L - u_L) / (S_L - S_star);
        double U_star_L_rho = coef;
        double U_star_L_u = coef * S_star;
        double U_star_L_v = coef * v_L;
        double U_star_L_E = coef * (E_L / rho_L + (S_star - u_L) * (S_star + P_L / (rho_L * (S_L - u_L))));
        
        F_star_rho = F_L_rho + S_L * (U_star_L_rho - rho_L);
        F_star_u = F_L_u + S_L * (U_star_L_u - rho_L * u_L);
        F_star_v = F_L_v + S_L * (U_star_L_v - rho_L * v_L);
        F_star_E = F_L_E + S_L * (U_star_L_E - E_L);
    } else {
        double coef = rho_R * (S_R - u_R) / (S_R - S_star);
        double U_star_R_rho = coef;
        double U_star_R_u = coef * S_star;
        double U_star_R_v = coef * v_R;
        double U_star_R_E = coef * (E_R / rho_R + (S_star - u_R) * (S_star + P_R / (rho_R * (S_R - u_R))));
        
        F_star_rho = F_R_rho + S_R * (U_star_R_rho - rho_R);
        F_star_u = F_R_u + S_R * (U_star_R_u - rho_R * u_R);
        F_star_v = F_R_v + S_R * (U_star_R_v - rho_R * v_R);
        F_star_E = F_R_E + S_R * (U_star_R_E - E_R);
    }

    F->rho = F_star_rho;
    if (dir == 0) { F->px = F_star_u; F->py = F_star_v; }
    else { F->py = F_star_u; F->px = F_star_v; }
    F->E = F_star_E;
}
