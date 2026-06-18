import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

# 1. Parameter parsing
parser = argparse.ArgumentParser(description="Compare SPH KH instability simulation data at different resolutions.")
parser.add_argument("--dirs", nargs='+', default=["sph_kh_2d_n32", "sph_kh_2d_n64", "sph_kh_2d_n128"], 
                    help="List of directories containing SPH output_*.csv files")
parser.add_argument("--labels", nargs='+', default=["SPH N=32", "SPH N=64", "SPH N=128"], 
                    help="List of labels corresponding to the SPH directories")
parser.add_argument("-o", "--output", type=str, default="kh_analysis_resolution.png", 
                    help="Output image filename (default: kh_analysis_resolution.png)")
args = parser.parse_args()

# Ensure dirs and labels match in length
if len(args.dirs) != len(args.labels):
    print("Warning: Number of SPH directories does not match number of labels. Using folder names as labels.")
    args.labels = [os.path.basename(d) for d in args.dirs]

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
        filename = os.path.basename(f)
        try:
            frame_idx = int(filename.split("_")[1].split(".")[0])
        except:
            frame_idx = idx
            
        t = frame_idx * dt_output
        times.append(t)
        
        try:
            df = pd.read_csv(f)
            m = df['m'] if 'm' in df.columns else df['mass']
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
        except Exception as e:
            print(f"Error reading file {f}: {e}")
            break
        
    return {
        "times": np.array(times),
        "Eky": np.array(Eky_list),
        "Ek": np.array(Ek_list),
        "U": np.array(U_list),
        "Etotal": np.array(Etotal_list),
        "Px": np.array(Px_list),
        "Py": np.array(Py_list)
    }

# Analyze SPH directories
data_runs = {}
for d, label in zip(args.dirs, args.labels):
    res = analyze_dir(d)
    if res is not None:
        data_runs[label] = res

if not data_runs:
    print("Error: No data could be processed. Exiting.")
    exit(1)

# Set up styling with larger fonts
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
})

fig, axs = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# Curated color palette
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

# --- Plot 1: KH Instability Growth Rate (ln(Eky)) ---
# Theoretical KH growth rate
rho1, rho2 = 1.0, 2.0
v1, v2 = -0.5, 0.5
L = 1.0
k = 2 * (2 * np.pi) / L
tau_KH = (rho1 + rho2) / (abs(v2 - v1) * k * np.sqrt(rho1 * rho2))
# Eky ~ vy^2, so Eky ~ exp(2t / tau_KH) -> ln(Eky) ~ 2t / tau_KH
expected_Eky_slope = 2.0 / tau_KH

ref_t = None
ref_Eky_log = None

color_idx = 0
for label, data in data_runs.items():
    color = colors[color_idx % len(colors)]
    axs[0].plot(data["times"], np.log(data["Eky"] + 1e-20), '-', linewidth=2.5, color=color, label=label)
    color_idx += 1
        
    # Get a reference point for the theoretical line (e.g., t=0.2)
    if ref_t is None and len(data["times"]) > 20:
        idx_ref = min(20, len(data["times"]) - 1)
        ref_t = data["times"][idx_ref]
        ref_Eky_log = np.log(data["Eky"][idx_ref] + 1e-20)

# Plot theoretical linear growth
if ref_t is not None:
    t_theory = np.linspace(0.0, 0.8, 100)
    # y - y1 = m(x - x1) -> y = m(x - x1) + y1
    Eky_theory_log = expected_Eky_slope * (t_theory - ref_t) + ref_Eky_log
    axs[0].plot(t_theory, Eky_theory_log, 'k:', linewidth=3.0, label="Analytic Linear Growth")

axs[0].set_title("Kelvin-Helmholtz Instability Growth Rate", fontweight='bold')
axs[0].set_ylabel("ln($E_{k,y}$)")
axs[0].grid(True, linestyle='--', alpha=0.6)
axs[0].legend(frameon=True, loc='lower right')
axs[0].set_ylim(bottom=ref_Eky_log - 4 if ref_Eky_log else -15) # zoom in reasonably

# --- Plot 2: Relative Energy Conservation Error (E(t) - E(0)) / E(0) ---
color_idx = 0
for label, data in data_runs.items():
    E0 = data["Etotal"][0]
    rel_error = (data["Etotal"] - E0) / E0
    color = colors[color_idx % len(colors)]
    axs[1].plot(data["times"], rel_error, '-', linewidth=2.5, color=color, label=label)
    color_idx += 1

axs[1].set_title("Relative Total Energy Error: $(E(t) - E(0)) / E(0)$", fontweight='bold')
axs[1].set_ylabel("Relative Error")
axs[1].grid(True, linestyle='--', alpha=0.6)
axs[1].legend(frameon=True)

# --- Plot 3: Specific Kinetic Energy (Ek) & Specific Internal Energy (U) ---
color_idx = 0
for label, data in data_runs.items():
    color = colors[color_idx % len(colors)]
    axs[2].plot(data["times"], data["Ek"], '-', linewidth=2.5, color=color, label=f"{label} Kinetic ($E_k$)")
    axs[2].plot(data["times"], data["U"], '--', linewidth=2.5, color=color, alpha=0.7, label=f"{label} Internal ($U$)")
    color_idx += 1

axs[2].set_title("Energy Partition: Kinetic vs. Internal Energy", fontweight='bold')
axs[2].set_xlabel("Time (t)")
axs[2].set_ylabel("Energy")
axs[2].grid(True, linestyle='--', alpha=0.6)
axs[2].legend(ncol=2, frameon=True)

plt.suptitle("Kelvin-Helmholtz Instability Resolution Comparison", fontsize=22, fontweight='bold', y=0.94)
plt.tight_layout(rect=[0, 0, 1, 0.93])

output_image = args.output
plt.savefig(output_image, dpi=200, bbox_inches='tight')
print(f"Comparison plot saved to {output_image}")
