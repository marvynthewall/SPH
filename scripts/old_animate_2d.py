import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os

parser = argparse.ArgumentParser(description="Create SPH animation from CSV files in a specific directory.")
parser.add_argument("-d", "--dir", type=str, required=True, help="Target directory containing output_*.csv files")
parser.add_argument("-o", "--outputfilename", type=str, default="sod_animation.gif", help="output gif filename")
args = parser.parse_args()

search_pattern = os.path.join(args.dir, "output_*.csv")

# Find all output csv files and sort them
files = sorted(glob.glob(search_pattern))
# files = sorted(glob.glob("output_*.csv"))
if not files:
    print("No output files found!")
    exit(1)

fig, ax = plt.subplots(figsize=(15, 3))

# Initialize the plot with the first frame
df = pd.read_csv(files[0])
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=5, vmin=0.0, vmax=1.2, cmap='viridis')
fig.colorbar(sc, ax=ax, label="Density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"SPH 2D Sod Shock Tube - {files[0]}")
ax.set_xlim(0, 5.0)
ax.set_ylim(0, 1.0)

# let the x, y axis in scale
ax.set_aspect('equal')

def update(frame):
    filename = files[frame]
    df = pd.read_csv(filename)
    
    # Update scatter plot data
    sc.set_offsets(df[["x", "y"]].values)
    sc.set_array(df["rho"].values)
    ax.set_title(f"SPH 2D Sod Shock Tube - {filename}")
    return sc,


def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

# Save as GIF
output_gif = args.outputfilename
print(f"Saving animation to {output_gif}...")
# ani.save(output_gif, writer='pillow', fps=10)
ani.save(output_gif, writer='pillow', fps=10, progress_callback=print_progress)
print("Done!")
