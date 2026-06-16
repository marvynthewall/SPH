import numpy as np
from scipy.optimize import fsolve
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os

def get_sod_exact(x_array, t, x0=7.5):
    """
    Calculate the exact analytical solution for the 1D Sod Shock Tube.
    Returns: rho_exact, u_exact, P_exact
    """
    gamma = 1.4
    g_m1 = gamma - 1.0
    g_p1 = gamma + 1.0
    
    rho_L, P_L, u_L = 1.0, 1.0, 0.0
    rho_R, P_R, u_R = 0.125, 0.1, 0.0
    
    if t == 0:
        rho_exact = np.where(x_array < x0, rho_L, rho_R)
        u_exact = np.where(x_array < x0, u_L, u_R)
        P_exact = np.where(x_array < x0, P_L, P_R)
        return rho_exact, u_exact, P_exact

    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)

    def shock_tube_equation(P):
        u_star_L = u_L + (2 * c_L / g_m1) * (1 - (P / P_L)**(g_m1 / (2 * gamma)))
        A = 2.0 / (g_p1 * rho_R)
        B = g_m1 / g_p1 * P_R
        u_star_R = u_R + (P - P_R) * np.sqrt(A / (P + B))
        return u_star_L - u_star_R

    P_star = fsolve(shock_tube_equation, 0.5)[0]
    u_star = u_L + (2 * c_L / g_m1) * (1 - (P_star / P_L)**(g_m1 / (2 * gamma)))
    rho_star_L = rho_L * (P_star / P_L)**(1 / gamma)
    rho_star_R = rho_R * ((P_star / P_R + g_m1 / g_p1) / (g_m1 / g_p1 * P_star / P_R + 1))
    c_star_L = np.sqrt(gamma * P_star / rho_star_L)

    v_shock = u_R + c_R * np.sqrt(g_p1 / (2 * gamma) * (P_star / P_R) + g_m1 / (2 * gamma))
    v_head = u_L - c_L
    v_tail = u_star - c_star_L
    v_contact = u_star

    rho_exact = np.zeros_like(x_array)
    u_exact = np.zeros_like(x_array)
    P_exact = np.zeros_like(x_array)
    
    for i, x in enumerate(x_array):
        xi = (x - x0) / t
        if xi < v_head:
            rho_exact[i] = rho_L
            u_exact[i] = u_L
            P_exact[i] = P_L
        elif v_head <= xi < v_tail:
            c = (2.0 / g_p1) * c_L + (g_m1 / g_p1) * (u_L - xi)
            rho_exact[i] = rho_L * (c / c_L)**(2 / g_m1)
            u_exact[i] = (2.0 / g_p1) * (c_L + g_m1 / 2.0 * u_L + xi)
            P_exact[i] = P_L * (c / c_L)**(2 * gamma / g_m1)
        elif v_tail <= xi < v_contact:
            rho_exact[i] = rho_star_L
            u_exact[i] = u_star
            P_exact[i] = P_star
        elif v_contact <= xi < v_shock:
            rho_exact[i] = rho_star_R
            u_exact[i] = u_star
            P_exact[i] = P_star
        else:
            rho_exact[i] = rho_R
            u_exact[i] = u_R
            P_exact[i] = P_R
            
    return rho_exact, u_exact, P_exact

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Sod Shock Tube Solutions")
    parser.add_argument("--sph1d_dir", type=str, default="sod_1d_m0.001_x15.0_t5.0")
    parser.add_argument("--sph2d_dir", type=str, default="sod_m0.001_x15.0_t5.0")
    parser.add_argument("--sph3d_dir", type=str, default="sod_3d_m0.001_x15.0_t5.0")
    parser.add_argument("--grid_dir", type=str, default="grid_sod_2d_output")
    parser.add_argument("-t", "--time", type=float, default=5.0)
    parser.add_argument("-x", "--xmax", type=float, default=15.0)
    parser.add_argument("-o", "--output", type=str, default="sod_comparison_all.png")
    args = parser.parse_args()

    file_idx = int(round(args.time / 0.01))
    
    def get_file(dir_path):
        return os.path.join(dir_path, f"output_{file_idx:04d}.csv")

    fig, axs = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    
    # 1. Exact Solution
    x_exact = np.linspace(0, args.xmax, 1000)
    rho_exact, u_exact, P_exact = get_sod_exact(x_exact, args.time, x0=args.xmax/2.0)
    
    axs[0].plot(x_exact, rho_exact, 'k-', lw=2, label="Exact")
    axs[1].plot(x_exact, u_exact, 'k-', lw=2, label="Exact")
    axs[2].plot(x_exact, P_exact, 'k-', lw=2, label="Exact")

    # Helper function to plot
    def plot_data(df, label, color, marker, ls='', alpha=1.0, s=2):
        # Calculate exact solution at the locations in df["x"]
        rho_ex, u_ex, P_ex = get_sod_exact(df["x"].values, args.time, x0=args.xmax/2.0)
        
        p_col = "P" if "P" in df.columns else "pressure"
        l1_rho = np.mean(np.abs(df["rho"].values - rho_ex))
        l1_u = np.mean(np.abs(df["vx"].values - u_ex))
        l1_P = np.mean(np.abs(df[p_col].values - P_ex))
        
        lbl_rho = f"{label} ($L_1$: {l1_rho:.4e})"
        lbl_u = f"{label} ($L_1$: {l1_u:.4e})"
        lbl_P = f"{label} ($L_1$: {l1_P:.4e})"

        if "Grid" in label:
            df = df.sort_values(by="x")
            axs[0].plot(df["x"], df["rho"], color=color, marker=marker, ls='-', ms=3, alpha=alpha, label=lbl_rho)
            axs[1].plot(df["x"], df["vx"], color=color, marker=marker, ls='-', ms=3, alpha=alpha, label=lbl_u)
            axs[2].plot(df["x"], df["P"], color=color, marker=marker, ls='-', ms=3, alpha=alpha, label=lbl_P)
        else:
            axs[0].scatter(df["x"], df["rho"], color=color, s=s, alpha=alpha, label=lbl_rho)
            axs[1].scatter(df["x"], df["vx"], color=color, s=s, alpha=alpha, label=lbl_u)
            if p_col in df.columns:
                axs[2].scatter(df["x"], df[p_col], color=color, s=s, alpha=alpha, label=lbl_P)

    # 2. Grid Solver
    grid_f = get_file(args.grid_dir)
    if os.path.exists(grid_f):
        df_grid = pd.read_csv(grid_f)
        plot_data(df_grid, "Grid (HLLC)", "blue", ".", alpha=0.5)
    
    # 3. 1D SPH
    sph1d_f = get_file(args.sph1d_dir)
    if os.path.exists(sph1d_f):
        df_1d = pd.read_csv(sph1d_f)
        plot_data(df_1d, "SPH 1D", "red", ".", alpha=0.5)

    # 4. 2D SPH
    sph2d_f = get_file(args.sph2d_dir)
    if os.path.exists(sph2d_f):
        df_2d = pd.read_csv(sph2d_f)
        # SPH 2D has many particles, plotting all might be slow/messy. 
        # We take a thin slice around y=0.5
        y_mid = df_2d["y"].max() / 2.0
        df_2d_slice = df_2d[(df_2d["y"] > y_mid - 0.05) & (df_2d["y"] < y_mid + 0.05)]
        if len(df_2d_slice) == 0: df_2d_slice = df_2d # fallback
        plot_data(df_2d_slice, "SPH 2D", "green", ".", s=3, alpha=0.3)

    # 5. 3D SPH
    sph3d_f = get_file(args.sph3d_dir)
    if os.path.exists(sph3d_f):
        df_3d = pd.read_csv(sph3d_f)
        y_mid = df_3d["y"].max() / 2.0
        z_mid = df_3d["z"].max() / 2.0
        df_3d_slice = df_3d[(df_3d["y"] > y_mid - 0.05) & (df_3d["y"] < y_mid + 0.05) & 
                            (df_3d["z"] > z_mid - 0.05) & (df_3d["z"] < z_mid + 0.05)]
        if len(df_3d_slice) == 0: df_3d_slice = df_3d # fallback
        plot_data(df_3d_slice, "SPH 3D", "orange", ".", alpha=0.1)

    axs[0].set_ylabel("Density (rho)", fontsize=12)
    axs[1].set_ylabel("Velocity (u)", fontsize=12)
    axs[2].set_ylabel("Pressure (P)", fontsize=12)
    axs[2].set_xlabel("Position (x)", fontsize=12)

    for ax in axs:
        ax.set_xlim(0, args.xmax)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.legend(loc='best', fontsize=9)

    plt.suptitle(f"Sod Shock Tube Comparison at t = {args.time:.2f}", fontsize=16)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Comparison plot saved to {args.output}")
