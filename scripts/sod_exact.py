# sod_exact.py
import numpy as np
from scipy.optimize import fsolve

def get_sod_exact_rho(x_array, t, x0=2.5):
    """
    計算標準 Sod Shock Tube 在時間 t 的解析解密度。
    假設 gamma = 1.4
    左側 (L): rho=1.0, P=1.0, u=0.0
    右側 (R): rho=0.125, P=0.1, u=0.0
    x0: 震波初始位置 (根據你的設定，如果長度是 5.0，通常初始在中線 2.5)
    """
    if t == 0:
        return np.where(x_array < x0, 1.0, 0.125)

    gamma = 1.4
    g_m1 = gamma - 1
    g_p1 = gamma + 1
    
    rho_L, P_L, u_L = 1.0, 1.0, 0.0
    rho_R, P_R, u_R = 0.125, 0.1, 0.0
    
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)

    # 用 fsolve 找出中間區域的壓力 P_star
    def shock_tube_equation(P):
        # 稀疏波 (Rarefaction wave) 關係式
        u_star_L = u_L + (2 * c_L / g_m1) * (1 - (P / P_L)**(g_m1 / (2 * gamma)))
        # 震波 (Shock wave) 關係式
        A = 2 / (g_p1 * rho_R)
        B = g_m1 / g_p1 * P_R
        u_star_R = u_R + (P - P_R) * np.sqrt(A / (P + B))
        return u_star_L - u_star_R

    P_star = fsolve(shock_tube_equation, 0.5)[0]
    u_star = u_L + (2 * c_L / g_m1) * (1 - (P_star / P_L)**(g_m1 / (2 * gamma)))
    rho_star_L = rho_L * (P_star / P_L)**(1 / gamma)
    rho_star_R = rho_R * ((P_star / P_R + g_m1 / g_p1) / (g_m1 / g_p1 * P_star / P_R + 1))
    c_star_L = np.sqrt(gamma * P_star / rho_star_L)

    # 計算波的傳遞速度
    v_shock = u_R + c_R * np.sqrt(g_p1 / (2 * gamma) * (P_star / P_R) + g_m1 / (2 * gamma))
    v_head = u_L - c_L
    v_tail = u_star - c_star_L
    v_contact = u_star

    rho_exact = np.zeros_like(x_array)
    for i, x in enumerate(x_array):
        xi = (x - x0) / t
        if xi < v_head:
            rho_exact[i] = rho_L
        elif v_head <= xi < v_tail:
            c = (2.0 / g_p1) * c_L + (g_m1 / g_p1) * (u_L - xi)
            rho_exact[i] = rho_L * (c / c_L)**(2 / g_m1)
        elif v_tail <= xi < v_contact:
            rho_exact[i] = rho_star_L
        elif v_contact <= xi < v_shock:
            rho_exact[i] = rho_star_R
        else:
            rho_exact[i] = rho_R
            
    return rho_exact
