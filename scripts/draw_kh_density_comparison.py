import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Define the directories and their labels
runs = {
    "32 x 32": "sph_kh_2d_n32",
    "64 x 64": "sph_kh_2d_n64",
    "128 x 128": "sph_kh_2d_n128",
    "256 x 256": "sph_kh_2d_n256"
}

fig, axs = plt.subplots(2, 2, figsize=(10, 10))
axs = axs.ravel()

# Try to find data and plot
for idx, (label, dir_name) in enumerate(runs.items()):
    ax = axs[idx]
    search_pattern = os.path.join(dir_name, "output_*.csv")
    files = sorted(glob.glob(search_pattern))
    
    if not files:
        ax.text(0.5, 0.5, f"No data found in\n{dir_name}", 
                ha='center', va='center', fontsize=12, color='gray')
        ax.set_title(label, fontsize=14, fontweight='bold')
        ax.axis('off')
        continue
        
    # Get the latest file
    latest_file = files[-1]
    filename = os.path.basename(latest_file)
    try:
        frame_idx = int(filename.split("_")[1].split(".")[0])
    except:
        frame_idx = len(files) - 1
        
    # Try to load time_log.csv if it exists
    time_log_path = os.path.join(dir_name, "time_log.csv")
    t = frame_idx * 0.01
    if os.path.exists(time_log_path):
        try:
            time_df = pd.read_csv(time_log_path)
            time_map = dict(zip(time_df["frame"], time_df["t"]))
            if frame_idx in time_map:
                t = time_map[frame_idx]
        except:
            pass
            
    print(f"Plotting {label} from {latest_file} (t = {t:.2f})")
    
    df = pd.read_csv(latest_file)
    
    # Scatter plot matching paper style
    # Marker size s can be adjusted based on resolution to make it look smooth
    if "32" in label:
        s_val = 25
    elif "64" in label:
        s_val = 8
    elif "128" in label:
        s_val = 3
    else:
        s_val = 1.0
        
    sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=s_val, zorder=2, vmin=0.9, vmax=2.1, cmap='jet')
    
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect('equal')
    
    # Add title with time info
    ax.set_title(f"{label} (t = {t:.2f})", fontsize=13, fontweight='bold')
    
    # Add a thin white border to subplots
    for spine in ax.spines.values():
        spine.set_color('white')
        spine.set_linewidth(1.5)
        
    # Keep it clean
    ax.set_xticks([])
    ax.set_yticks([])

# Add a single colorbar for all plots
fig.subplots_adjust(right=0.85, wspace=0.15, hspace=0.15)
cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
fig.colorbar(sc, cax=cbar_ax, label="Density (rho)")

plt.suptitle("SPH Kelvin-Helmholtz Instability Resolution Study", fontsize=16, fontweight='bold', y=0.95)

output_filename = "kh_resolution_comparison.png"
plt.savefig(output_filename, dpi=200, bbox_inches='tight')
print(f"Comparison plot saved to {output_filename}")
