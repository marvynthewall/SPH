import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

# 1. Parameter parsing
parser = argparse.ArgumentParser(description="Compare SPH and Grid KH instability simulation data.")
parser.add_argument("--sph_dir", type=str, default="sph_kh_2d_output", 
                    help="Directory containing SPH output_*.csv files")
parser.add_argument("--grid_dir", type=str, default="grid_kh_2d_output", 
                    help="Directory containing Grid output_*.csv files")
parser.add_argument("-o", "--output", type=str, default="kh_comparison.png", 
                    help="Output image filename (default: kh_comparison.png)")
args = parser.parse_args()

def analyze_dir(dir_path):
    search_pattern = os.path.join(dir_path, "output_*.csv")
    files = sorted(glob.glob(search_pattern))
    
    if not files:
        print(f"Warning: No files found matching '{search_pattern}' in '{dir_path}'")
        return None

    print(f"Processing {len(files)} files in '{dir_path}'...")
    times = []
    Eky_list = []
    Ek_list = []
    U_list = []
    Etotal_list = []
    Px_list = []
    Py_list = []
    
    dt_output = 0.01  # Output interval
    
    for idx, f in enumerate(files):
        df = pd.read_csv(f)
        t = idx * dt_output
        times.append(t)
        
        m = df['m']
        vx = df['vx']
        vy = df['vy']
        u = df['u']
        
        Eky = np.sum(0.5 * m * vy**2)
        Eky_list.append(Eky)
        
        Ekx = np.sum(0.5 * m * vx**2)
        Ek = Ekx + Eky
        Ek_list.append(Ek)
        
        U = np.sum(m * u)
        U_list.append(U)
        
        Etotal = Ek + U
        Etotal_list.append(Etotal)
        
        Px = np.sum(m * vx)
        Py = np.sum(m * vy)
        Px_list.append(Px)
        Py_list.append(Py)
        
    return {
        "times": np.array(times),
        "Eky": np.array(Eky_list),
        "Ek": np.array(Ek_list),
        "U": np.array(U_list),
        "Etotal": np.array(Etotal_list),
        "Px": np.array(Px_list),
        "Py": np.array(Py_list)
    }

print("Starting analysis...")
sph_data = analyze_dir(args.sph_dir)
grid_data = analyze_dir(args.grid_dir)

if sph_data is None and grid_data is None:
    print("Error: No data found for both SPH and Grid. Exiting.")
    exit(1)

fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# --- Plot 1: KH Instability Growth Rate (ln(Eky)) ---
if sph_data is not None:
    axs[0].plot(sph_data["times"], np.log(sph_data["Eky"] + 1e-20), 'r-', linewidth=2, label="SPH 2D")
if grid_data is not None:
    axs[0].plot(grid_data["times"], np.log(grid_data["Eky"] + 1e-20), 'b--', linewidth=2, label="Grid (MUSCL)")
axs[0].set_title("Kelvin-Helmholtz Instability Growth Rate\n(Linear phase shows constant slope)", fontsize=14)
axs[0].set_xlabel("Time (t)")
axs[0].set_ylabel("ln(E_{k,y})")
axs[0].grid(True, linestyle='--', alpha=0.7)
axs[0].legend()

# --- Plot 2: Energy Conservation Check ---
if grid_data is not None:
    axs[1].plot(grid_data["times"], grid_data["Ek"], 'b-', alpha=0.5, label='Grid Kinetic ($E_k$)')
    axs[1].plot(grid_data["times"], grid_data["U"], 'r-', alpha=0.5, label='Grid Internal ($U$)')
    axs[1].plot(grid_data["times"], grid_data["Etotal"], 'k-', linewidth=2, label='Grid Total ($E_{total}$)')
if sph_data is not None:
    axs[1].plot(sph_data["times"], sph_data["Ek"], 'b--', alpha=0.8, label='SPH Kinetic ($E_k$)')
    axs[1].plot(sph_data["times"], sph_data["U"], 'r--', alpha=0.8, label='SPH Internal ($U$)')
    axs[1].plot(sph_data["times"], sph_data["Etotal"], 'k--', linewidth=2, label='SPH Total ($E_{total}$)')
axs[1].set_title("Energy Conservation Check", fontsize=14)
axs[1].set_xlabel("Time (t)")
axs[1].set_ylabel("Energy")
axs[1].legend(ncol=2)
axs[1].grid(True, linestyle='--', alpha=0.7)

# --- Plot 3: Momentum Conservation Check ---
if grid_data is not None:
    axs[2].plot(grid_data["times"], grid_data["Px"], 'm-', alpha=0.6, label='Grid Total $P_x$')
    axs[2].plot(grid_data["times"], grid_data["Py"], 'c-', alpha=0.6, label='Grid Total $P_y$')
if sph_data is not None:
    axs[2].plot(sph_data["times"], sph_data["Px"], 'm--', alpha=0.9, label='SPH Total $P_x$')
    axs[2].plot(sph_data["times"], sph_data["Py"], 'c--', alpha=0.9, label='SPH Total $P_y$')
axs[2].set_title("Momentum Conservation Check", fontsize=14)
axs[2].set_xlabel("Time (t)")
axs[2].set_ylabel("Momentum")
axs[2].legend(ncol=2)
axs[2].grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
output_image = args.output
plt.savefig(output_image, dpi=150)
print(f"Comparison plot saved to {output_image}")
