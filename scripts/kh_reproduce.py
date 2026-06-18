import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

# ==========================================
# 1. Parameter Parsing
# ==========================================
parser = argparse.ArgumentParser(description="Reproduce Springel (2011) Fig 8: Mode Amplitude & Velocity Dispersion for SPH KHI.")
parser.add_argument("--dirs", nargs='+', default=["sph_kh_2d_n32", "sph_kh_2d_n64", "sph_kh_2d_n128"], 
                    help="List of directories containing SPH output_*.csv files")
parser.add_argument("--labels", nargs='+', default=["32x32", "64x64", "128x128"], 
                    help="List of labels corresponding to the SPH directories")
parser.add_argument("-o", "--output", type=str, default="kh_springel_fig8.png", 
                    help="Output image filename (default: kh_springel_fig8.png)")
args = parser.parse_args()

if len(args.dirs) != len(args.labels):
    print("Warning: Number of directories does not match labels. Using folder names.")
    args.labels = [os.path.basename(d) for d in args.dirs]

# Physical parameters from Springel 2011
rho1, rho2 = 1.0, 2.0
v1, v2 = -0.5, 0.5
L = 1.0
MODE = 2.0  # Wavenumber k = 2 * (2pi/L)
k = MODE * (2.0 * np.pi) / L
tau_KH = (rho1 + rho2) / (abs(v2 - v1) * k * np.sqrt(rho1 * rho2))

def analyze_dir(dir_path):
    search_pattern = os.path.join(dir_path, "output_*.csv")
    files = sorted(glob.glob(search_pattern))
    
    if not files:
        print(f"Warning: No files found matching '{search_pattern}' in '{dir_path}'")
        return None

    print(f"Processing {len(files)} files in '{dir_path}'...")
    times = []
    mode_amp_list = []
    sigma_vx_list = []
    
    dt_output = 0.01  # Assuming output interval is 0.01
    
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
            # Ensure coordinates and velocities exist
            x = df['x'].values
            y = df['y'].values
            vx = df['vx'].values
            vy = df['vy'].values
            m = df['m'].values if 'm' in df.columns else df['mass'].values
            
            # -------------------------------------------------------------
            # 1. Calculate Mode Amplitude using Direct Fourier Projection
            # This is mathematically equivalent to picking the amplitude 
            # of wavenumber k from an FFT, but perfect for scattered data.
            # -------------------------------------------------------------
            vy_cos = np.average(vy * np.cos(k * x), weights=m)
            vy_sin = np.average(vy * np.sin(k * x), weights=m)
            # Amplitude of sine wave is 2 * magnitude of the Fourier component
            mode_amp = 2.0 * np.sqrt(vy_cos**2 + vy_sin**2)
            mode_amp_list.append(mode_amp)
            
            # -------------------------------------------------------------
            # 2. Calculate Velocity Dispersion at the Midplane
            # Springel: "velocity dispersion in x-direction for a thin layer
            # close to the midplane of the box"
            # -------------------------------------------------------------
            # Midplane is at y = 0.5. We take a thin layer of +/- 0.025
            midplane_mask = np.abs(y - 0.5) < 0.025
            if np.sum(midplane_mask) > 5: # Ensure we have enough particles
                vx_layer = vx[midplane_mask]
                sigma_vx = np.std(vx_layer)
            else:
                sigma_vx = np.nan
            sigma_vx_list.append(sigma_vx)
            
        except Exception as e:
            print(f"Error reading file {f}: {e}")
            break
        
    return {
        "times": np.array(times),
        "mode_amp": np.array(mode_amp_list),
        "sigma_vx": np.array(sigma_vx_list)
    }

# Process all directories
data_runs = {}
for d, label in zip(args.dirs, args.labels):
    res = analyze_dir(d)
    if res is not None:
        data_runs[label] = res

if not data_runs:
    print("Error: No data could be processed. Exiting.")
    exit(1)

# ==========================================
# 2. Plotting (Reproducing Paper's Style)
# ==========================================
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 13,
})

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Paper's line colors usually match specific resolutions: 
# Red(32), Blue(64), Magenta/Purple(128), Cyan(256), Green(512)
colors = ['#d62728', '#1f77b4', '#9467bd', '#17becf', '#2ca02c']

# ---------------------------------------------
# Left Panel: Mode Amplitude
# ---------------------------------------------
color_idx = 0
for label, data in data_runs.items():
    color = colors[color_idx % len(colors)]
    axs[0].plot(data["times"], data["mode_amp"], '-', linewidth=1.5, color=color, label=label)
    color_idx += 1

# Plot Analytic Linear Growth
# v_y(t) = v_0 * exp(t / tau_KH)
# Initial amplitude is v0 = 0.01
t_theory = np.linspace(0.0, 1.0, 100)
amp_theory = 0.01 * np.exp(t_theory / tau_KH)
axs[0].plot(t_theory, amp_theory, 'k--', linewidth=1.5, label='Analytic linear growth')

axs[0].set_yscale('log')
axs[0].set_xlabel('$t$')
axs[0].set_ylabel('mode amplitude')
axs[0].set_xlim(0.0, 1.0)
axs[0].set_ylim(0.001, 0.2)
axs[0].grid(True, which="both", ls="--", alpha=0.5)

# ---------------------------------------------
# Right Panel: Velocity Dispersion (Noise)
# ---------------------------------------------
color_idx = 0
for label, data in data_runs.items():
    color = colors[color_idx % len(colors)]
    axs[1].plot(data["times"], data["sigma_vx"], '-', linewidth=1.5, color=color, label=label)
    color_idx += 1

axs[1].set_yscale('log')
axs[1].set_xlabel('$t$')
axs[1].set_ylabel(r'$\sigma_{V_x}$')
axs[1].set_xlim(0.0, 1.0)
axs[1].set_ylim(0.001, 0.2)
axs[1].grid(True, which="both", ls="--", alpha=0.5)

# Adjust legends
axs[0].legend(loc='upper left', frameon=False)
axs[1].legend(loc='upper left', frameon=False)

plt.tight_layout()
output_image = args.output
plt.savefig(output_image, dpi=200, bbox_inches='tight')
print(f"Successfully generated {output_image}")