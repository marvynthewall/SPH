import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import time

# Initialize argument parser
parser = argparse.ArgumentParser(description="Plot SPH simulation density distribution")

# Add arguments
parser.add_argument("-f", "--filename", type=str, default="output_0000.csv",
                    help="CSV filename to read (default: output_0000.csv)")
parser.add_argument("-x", "--x", type=float, default=15.0,
                    help="size of x direction (default: 15.0)")

# 新增參數：決定是否要畫 h circles (預設為不畫，加上 -c 或 --circle 才會畫)
parser.add_argument("-c", "--circle", action="store_true", 
                    help="Draw smoothing length (h) circles (Warning: may be slow for large datasets)")

# Parse arguments
args = parser.parse_args()

# Read data using the parsed filename
print(f"[{time.strftime('%H:%M:%S')}] Reading file: {args.filename} ...")
df = pd.read_csv(args.filename)
print(f"[{time.strftime('%H:%M:%S')}] Loaded {len(df)} particles.")

# Create figure and axis
fig, ax = plt.subplots(figsize=(20,6))

# ---------------------------------------------------------
# 1. Draw semi-transparent circles for smoothing length (h)
# ---------------------------------------------------------
if args.circle:
    print(f"[{time.strftime('%H:%M:%S')}] Generating h circles...")
    circles = [Circle((x, y), r) for x, y, r in zip(df['x'], df['y'], df['h'])]
    
    print(f"[{time.strftime('%H:%M:%S')}] Adding PatchCollection to axis...")
    patch_col = PatchCollection(circles, alpha=0.15, facecolor='gray', edgecolor='none', zorder=1)
    ax.add_collection(patch_col)
else:
    print(f"[{time.strftime('%H:%M:%S')}] Skipping h circles (use --circle to enable).")

# ---------------------------------------------------------
# 2. Draw particle center points
# ---------------------------------------------------------
print(f"[{time.strftime('%H:%M:%S')}] Setting up scatter plot...")
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=0.2, zorder=2, cmap='viridis')

# Add colorbar and labels
fig.colorbar(sc, ax=ax, label="density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")

# Set the limits of x and y axes
ax.set_xlim(0.0, args.x)
ax.set_ylim(0.0, 1.0)

# Set title dynamically based on whether circles are drawn
title_str = f"SPH Density & Smoothing Length: {args.filename}" if args.circle else f"SPH Density: {args.filename}"
plt.title(title_str)

output_png = args.filename.replace(".csv", ".png")

# 提示使用者正在進行最耗時的渲染步驟
print(f"[{time.strftime('%H:%M:%S')}] Rendering and saving plot to {output_png}...")
start_render = time.time()

plt.savefig(output_png, bbox_inches='tight', dpi=150)

render_time = time.time() - start_render
print(f"[{time.strftime('%H:%M:%S')}] Done! Rendering took {render_time:.2f} seconds.")

# 最後明確印出輸出的檔名
print(f"=========================================")
print(f"Output file successfully saved as: {output_png}")
print(f"=========================================")
