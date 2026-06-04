import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(
    description="Create 2D x-y field animation from 3D SPH output_*.csv files."
)

parser.add_argument("-d", "--dir", type=str, required=True,
                    help="Target directory containing output_*.csv files")

parser.add_argument("-o", "--outputfilename", type=str,
                    default="sod3d_2d_fields",
                    help="Output filename")

parser.add_argument("-x", "--x_lim", type=float, default=5.0,
                    help="x direction limit")

parser.add_argument("-y", "--y_lim", type=float, default=1.0,
                    help="y direction limit")

parser.add_argument("-z", "--z_lim", type=float, default=1.0,
                    help="z direction limit")

parser.add_argument("--z_slice", type=float, default=0.5,
                    help="z position of the 2D slice")

parser.add_argument("--slice_width", type=float, default=0.05,
                    help="Thickness of the z slice")

parser.add_argument("-f", "--format", type=str, default="mp4",
                    help="File format: mp4 or gif")

parser.add_argument("--fps", type=int, default=10,
                    help="Frames per second")

parser.add_argument("--point_size", type=float, default=8.0,
                    help="Particle marker size")

parser.add_argument("--gamma", type=float, default=1.4,
                    help="Gas gamma, used if pressure column is missing")

args = parser.parse_args()

# ------------------------------------------------------------
# 1. Find output files
# ------------------------------------------------------------
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}")
    exit(1)

print(f"Found {len(files)} files.")

# ------------------------------------------------------------
# 2. Read time log if available
# ------------------------------------------------------------
time_log_path = os.path.join(args.dir, "time_log.csv")
time_map = {}

if os.path.exists(time_log_path):
    time_df = pd.read_csv(time_log_path)
    if "frame" in time_df.columns and "t" in time_df.columns:
        for _, row in time_df.iterrows():
            time_map[int(row["frame"])] = float(row["t"])

# ------------------------------------------------------------
# 3. Helper: load one frame and extract z slice
# ------------------------------------------------------------
def load_slice(filepath):
    df = pd.read_csv(filepath)

    required_cols = ["x", "y", "z", "rho", "vx"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"CSV file missing required column: {col}")

    # pressure
    if "pressure" in df.columns:
        pressure = df["pressure"].to_numpy()
    elif "u" in df.columns:
        pressure = (args.gamma - 1.0) * df["rho"].to_numpy() * df["u"].to_numpy()
    else:
        raise ValueError("CSV needs either 'pressure' column or 'u' column.")

    df = df.copy()
    df["pressure_plot"] = pressure

    # z slice
    z_min = args.z_slice - 0.5 * args.slice_width
    z_max = args.z_slice + 0.5 * args.slice_width

    df_slice = df[(df["z"] >= z_min) & (df["z"] <= z_max)]

    if len(df_slice) == 0:
        raise ValueError(
            f"No particles found in z slice: z={args.z_slice}, width={args.slice_width}. "
            f"Try increasing --slice_width."
        )

    return df_slice

# ------------------------------------------------------------
# 4. Initialize first frame
# ------------------------------------------------------------
df0 = load_slice(files[0])

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

common_cmap = "viridis"

sc_rho = axes[0].scatter(
    df0["x"], df0["y"],
    c=df0["rho"],
    s=args.point_size,
    cmap=common_cmap,
    vmin=0.0,
    vmax=1.2
)

sc_vx = axes[1].scatter(
    df0["x"], df0["y"],
    c=df0["vx"],
    s=args.point_size,
    cmap=common_cmap,
    vmin=-0.5,
    vmax=1.5
)

sc_p = axes[2].scatter(
    df0["x"], df0["y"],
    c=df0["pressure_plot"],
    s=args.point_size,
    cmap=common_cmap,
    vmin=0.0,
    vmax=1.2
)

cbar0 = fig.colorbar(sc_rho, ax=axes[0])
cbar0.set_label(r"$\rho$")

cbar1 = fig.colorbar(sc_vx, ax=axes[1])
cbar1.set_label(r"$v_x$")

cbar2 = fig.colorbar(sc_p, ax=axes[2])
cbar2.set_label(r"$P$")

axes[0].set_ylabel("y")
axes[1].set_ylabel("y")
axes[2].set_ylabel("y")
axes[2].set_xlabel("x")

axes[0].set_title("Density")
axes[1].set_title("x-velocity")
axes[2].set_title("Pressure")

for ax in axes:
    ax.set_xlim(0.0, args.x_lim)
    ax.set_ylim(0.0, args.y_lim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)

main_title = fig.suptitle(
    f"SPH 3D Sod Shock Tube 2D Slice | z={args.z_slice:.3f}"
)

plt.tight_layout()

# ------------------------------------------------------------
# 5. Update animation
# ------------------------------------------------------------
def update(frame):
    filepath = files[frame]
    df = load_slice(filepath)

    xy = np.column_stack([df["x"].to_numpy(), df["y"].to_numpy()])

    sc_rho.set_offsets(xy)
    sc_rho.set_array(df["rho"].to_numpy())

    sc_vx.set_offsets(xy)
    sc_vx.set_array(df["vx"].to_numpy())

    sc_p.set_offsets(xy)
    sc_p.set_array(df["pressure_plot"].to_numpy())

    if frame in time_map:
        main_title.set_text(
            f"SPH 3D Sod Shock Tube 2D Slice | z={args.z_slice:.3f}, t={time_map[frame]:.4f}"
        )
    else:
        current_filename = os.path.basename(filepath)
        main_title.set_text(
            f"SPH 3D Sod Shock Tube 2D Slice | z={args.z_slice:.3f}, {current_filename}"
        )

    return sc_rho, sc_vx, sc_p, main_title

# ------------------------------------------------------------
# 6. Progress display
# ------------------------------------------------------------
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# ------------------------------------------------------------
# 7. Create animation
# ------------------------------------------------------------
ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=100,
    blit=False
)

output_file = args.outputfilename + "." + args.format

print(f"Saving animation to {output_file} ...")

if args.format == "gif":
    ani.save(
        output_file,
        writer="pillow",
        fps=args.fps,
        progress_callback=print_progress
    )

elif args.format == "mp4":
    ani.save(
        output_file,
        writer="ffmpeg",
        fps=args.fps,
        progress_callback=print_progress
    )

else:
    print("Wrong file format. Please use mp4 or gif.")
    exit(1)

print("\nDone!")


# python animate_sod3d_2d_fields.py \
#   -d ../sod3d_m0.001_x5.0_y1.0_z1.0_t5.0 \
#   -o sod3d_2d_fields \
#   -x 5.0 \
#   -y 1.0 \
#   -z 1.0 \
#   --z_slice 0.5 \
#   --slice_width 0.05 \
#   -f mp4 \
#   --fps 5