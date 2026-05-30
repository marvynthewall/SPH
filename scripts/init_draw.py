import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

# Initialize argument parser
parser = argparse.ArgumentParser(description="Plot SPH simulation density distribution")

# Add -f argument with a default value
parser.add_argument("-f", "--filename", type=str, default="output_0000.csv", 
                    help="CSV filename to read (default: output_0000.csv)")

# Parse arguments
args = parser.parse_args()

# Read data using the parsed filename
print(f"Reading file: {args.filename} ...")
df = pd.read_csv(args.filename)

# Create figure and axis
fig, ax = plt.subplots()

# ---------------------------------------------------------
# 1. Draw semi-transparent circles for smoothing length (h)
# ---------------------------------------------------------
# Use zip to quickly combine x, y, and h into a list of Circle objects
circles = [Circle((x, y), r) for x, y, r in zip(df['x'], df['y'], df['h'])]

# Pack all circles into a PatchCollection for better performance.
# alpha: Transparency (0.0 to 1.0)
# facecolor: Fill color
# edgecolor: Border color ('none' means no border)
# zorder: Layer order. Lower number means it's placed behind. Set to 1 for background.
patch_col = PatchCollection(circles, alpha=0.15, facecolor='gray', edgecolor='none', zorder=1)
ax.add_collection(patch_col)


# ---------------------------------------------------------
# 2. Draw particle center points
# ---------------------------------------------------------
# Set zorder=2 to ensure points are drawn on top of the circles (zorder=1)
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=10, zorder=2, cmap='viridis')

# Add colorbar and labels
fig.colorbar(sc, ax=ax, label="density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")

# Set the limits of x and y axes to 0.0 ~ 1.0
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)

# Set title with the filename
plt.title(f"SPH Density & Smoothing Length: {args.filename}") 

plt.show()
