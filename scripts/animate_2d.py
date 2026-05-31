import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob

# Find all output csv files and sort them
files = sorted(glob.glob("output_*.csv"))
if not files:
    print("No output files found!")
    exit(1)

fig, ax = plt.subplots(figsize=(8, 4))

# Initialize the plot with the first frame
df = pd.read_csv(files[0])
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=5, vmin=0.0, vmax=1.2, cmap='viridis')
fig.colorbar(sc, ax=ax, label="Density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"SPH 2D Sod Shock Tube - {files[0]}")
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.0)

def update(frame):
    filename = files[frame]
    df = pd.read_csv(filename)
    
    # Update scatter plot data
    sc.set_offsets(df[["x", "y"]].values)
    sc.set_array(df["rho"].values)
    ax.set_title(f"SPH 2D Sod Shock Tube - {filename}")
    return sc,

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

# Save as GIF
output_gif = "sod_animation.gif"
print(f"Saving animation to {output_gif}...")
ani.save(output_gif, writer='pillow', fps=10)
print("Done!")
