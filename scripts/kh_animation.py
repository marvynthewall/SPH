import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os

# 1. Parameter parsing
parser = argparse.ArgumentParser(description="Create KH instability animation.")
parser.add_argument("-d", "--dir", type=str, default=".", 
                    help="Target directory containing output_*.csv files (default: current directory)")
parser.add_argument("-o", "--outputfilename", type=str, default="kh_animation", 
                    help="Output filename without extension (default: kh_animation)")
parser.add_argument("-f", "--format", type=str, default="gif", 
                    help="File format, mp4 or gif (default: gif)")
args = parser.parse_args()

# 2. Find and sort CSV files
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}!")
    exit(1)

print(f"Found {len(files)} files. Initializing plot...")

# 3. Create figure and initial setup
fig, ax = plt.subplots(figsize=(8, 8))

# Read the first frame
df = pd.read_csv(files[0])

# --- Plot initial scatter ---
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=5, zorder=2, vmin=0.5, vmax=2.5, cmap='jet')

# Formatting
fig.colorbar(sc, ax=ax, label="Density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.0)

# Equal aspect ratio
ax.set_aspect('equal')

first_filename = os.path.basename(files[0])
title = ax.set_title(f"SPH KH Instability - {first_filename}")

# 4. Animation update function
def update(frame):
    filepath = files[frame]
    df = pd.read_csv(filepath)
    
    # Update scatter plot
    sc.set_offsets(df[["x", "y"]].values)
    sc.set_array(df["rho"].values)
    
    # Update title
    current_filename = os.path.basename(filepath)
    title.set_text(f"SPH KH Instability - {current_filename}")
    
    return sc, title

# 5. Progress bar callback
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 6. Generate and save animation
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

if args.format == 'gif':
    output_file = args.outputfilename + ".gif"
    print(f"Saving animation to {output_file} ...")
    ani.save(output_file, writer='pillow', fps=10, progress_callback=print_progress)
elif args.format == 'mp4':
    output_file = args.outputfilename + ".mp4"
    print(f"Saving animation to {output_file} ...")
    ani.save(output_file, writer='ffmpeg', fps=10, progress_callback=print_progress)
else:
    print(f"Error: Unsupported format '{args.format}'. Please use 'gif' or 'mp4'.")

print("\nDone!")
