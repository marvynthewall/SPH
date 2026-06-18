import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# Import the exact solver
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from sod_exact import get_sod_exact

def get_last_csv(dir_path):
    search_pattern = os.path.join(dir_path, "output_*.csv")
    files = sorted(glob.glob(search_pattern))
    if not files:
        return None
    return files[-1] # Return the last file which should be t_end

def calculate_l1_error(csv_path, target_time, xmax):
    df = pd.read_csv(csv_path)
    # Extract data near midplane for 2D/3D to avoid boundary artifacts
    if 'y' in df.columns:
        y_mid = df["y"].max() / 2.0
        mask = (df["y"] > y_mid - 0.05) & (df["y"] < y_mid + 0.05)
        if 'z' in df.columns:
            z_mid = df["z"].max() / 2.0
            mask = mask & (df["z"] > z_mid - 0.05) & (df["z"] < z_mid + 0.05)
        df_slice = df[mask]
        if len(df_slice) == 0: df_slice = df
    else:
        df_slice = df

    rho_sim = df_slice["rho"].values
    x_pos = df_slice["x"].values

    # Calculate exact solution at these specific particle/cell coordinates
    rho_exact, _, _ = get_sod_exact(x_pos, target_time, x0=xmax/2.0)
    
    # L1 Error
    l1_error = np.mean(np.abs(rho_sim - rho_exact))
    return l1_error

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=str, default="convergence_data")
    parser.add_argument("-t", "--time", type=float, default=1.0)
    parser.add_argument("-x", "--xmax", type=float, default=15.0)
    parser.add_argument("-o", "--output", type=str, default="sod_convergence.png")
    args = parser.parse_args()

    results = {
        "Grid 1D": {"x": [], "y": []},
        "SPH 1D": {"x": [], "y": []},
        "SPH 2D": {"x": [], "y": []},
        "SPH 3D": {"x": [], "y": []}
    }

    print("=== Calculating L1 Errors ===")
    
    # 1. Grid 1D
    grid_dirs = sorted(glob.glob(os.path.join(args.datadir, "grid_n*")))
    for d in grid_dirs:
        try:
            nx = int(os.path.basename(d).replace("grid_n", ""))
            grid_csv = get_last_csv(d)
            m_grid = args.xmax / nx  # Equivalent 1D mass (since L/nx = dx, and dx = m/rho)
            if grid_csv:
                err = calculate_l1_error(grid_csv, args.time, args.xmax)
                results["Grid 1D"]["x"].append(m_grid)
                results["Grid 1D"]["y"].append(err)
                print(f"Grid 1D (nx={nx}, equiv m={m_grid:.5f}): L1 = {err:.5e}")
        except ValueError:
            continue

    # 2. SPH 1D
    sph1d_dirs = sorted(glob.glob(os.path.join(args.datadir, "sph1d_m*")))
    for d in sph1d_dirs:
        try:
            m = float(os.path.basename(d).replace("sph1d_m", ""))
            sph1d_csv = get_last_csv(d)
            if sph1d_csv:
                err = calculate_l1_error(sph1d_csv, args.time, args.xmax)
                results["SPH 1D"]["x"].append(m)
                results["SPH 1D"]["y"].append(err)
                print(f"SPH 1D (m={m}): L1 = {err:.5e}")
        except ValueError:
            continue

    # 3. SPH 2D
    sph2d_dirs = sorted(glob.glob(os.path.join(args.datadir, "sph2d_m*")))
    for d in sph2d_dirs:
        try:
            m = float(os.path.basename(d).replace("sph2d_m", ""))
            sph2d_csv = get_last_csv(d)
            if sph2d_csv:
                err = calculate_l1_error(sph2d_csv, args.time, args.xmax)
                results["SPH 2D"]["x"].append(m)
                results["SPH 2D"]["y"].append(err)
                print(f"SPH 2D (m={m}): L1 = {err:.5e}")
        except ValueError:
            continue

    # 4. SPH 3D
    sph3d_dirs = sorted(glob.glob(os.path.join(args.datadir, "sph3d_m*")))
    for d in sph3d_dirs:
        try:
            m = float(os.path.basename(d).replace("sph3d_m", ""))
            sph3d_csv = get_last_csv(d)
            if sph3d_csv:
                err = calculate_l1_error(sph3d_csv, args.time, args.xmax)
                results["SPH 3D"]["x"].append(m)
                results["SPH 3D"]["y"].append(err)
                print(f"SPH 3D (m={m}): L1 = {err:.5e}")
        except ValueError:
            continue

    # Plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {"Grid 1D": "blue", "SPH 1D": "red", "SPH 2D": "green", "SPH 3D": "orange"}
    markers = {"Grid 1D": "s", "SPH 1D": "o", "SPH 2D": "^", "SPH 3D": "D"}

    for name, data in results.items():
        if len(data["x"]) > 1:
            label = name
            if name == "Grid 1D":
                label = "Grid 1D (aligned grid number with 1D SPH particle number)"
            ax.loglog(data["x"], data["y"], marker=markers[name], color=colors[name], 
                      linewidth=2, markersize=8, label=label)

    # Reference slopes (y = C * m^p)
    if len(results["Grid 1D"]["x"]) > 1:
        x_ref = np.array([results["Grid 1D"]["x"][0], results["Grid 1D"]["x"][-1]])
        y_ref_1 = results["Grid 1D"]["y"][0] * (x_ref / x_ref[0])**1.0
        ax.loglog(x_ref, y_ref_1, 'k--', alpha=0.5, label=r"1st Order $\mathcal{O}(m)$")

    ax.set_xlabel('Particle Mass ($m$)', fontsize=14)
    ax.set_ylabel('$L_1$ Error of Density ($\\rho$)', fontsize=14)
    ax.set_title('Sod Shock Tube Convergence Test', fontsize=16)
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=12)

    # Invert x-axis so higher resolution (smaller m) is on the right
    ax.invert_xaxis()

    plt.tight_layout()
    plt.savefig(args.output, dpi=200)
    print(f"\nConvergence plot saved to {args.output}")

if __name__ == "__main__":
    main()
