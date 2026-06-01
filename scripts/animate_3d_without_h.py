import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os
import numpy as np

# 1. 參數解析
parser = argparse.ArgumentParser(
    description="Create 3D SPH animation from output_*.csv files."
)

parser.add_argument("-d", "--dir", type=str, required=True,
                    help="Target directory containing output_*.csv files")

parser.add_argument("-o", "--outputfilename", type=str,
                    default="sod3d_animation",
                    help="Output filename")

parser.add_argument("-x", "--x_lim", type=float, default=5.0,
                    help="x direction limit")

parser.add_argument("-y", "--y_lim", type=float, default=1.0,
                    help="y direction limit")

parser.add_argument("-z", "--z_lim", type=float, default=1.0,
                    help="z direction limit")

parser.add_argument("-f", "--format", type=str, default="mp4",
                    help="File format: mp4 or gif")

parser.add_argument("--fps", type=int, default=10,
                    help="Frames per second")

parser.add_argument("--point_size", type=float, default=8.0,
                    help="Particle marker size")

parser.add_argument("--rho_min", type=float, default=0.0,
                    help="Minimum density for color scale")

parser.add_argument("--rho_max", type=float, default=1.2,
                    help="Maximum density for color scale")

args = parser.parse_args()

# 2. 找 output_*.csv
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}")
    exit(1)

print(f"Found {len(files)} files. Initializing 3D plot...")

# 3. 讀第一幀
df = pd.read_csv(files[0])

required_cols = ["x", "y", "z", "rho"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"CSV file missing required column: {col}")

# 4. 建立 3D figure
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection="3d")

sc = ax.scatter(
    df["x"],
    df["y"],
    df["z"],
    c=df["rho"],
    s=args.point_size,
    cmap="viridis",
    vmin=args.rho_min,
    vmax=args.rho_max,
    depthshade=True
)

cbar = fig.colorbar(sc, ax=ax, pad=0.1)
cbar.set_label("Density")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.set_xlim(0, args.x_lim)
ax.set_ylim(0, args.y_lim)
ax.set_zlim(0, args.z_lim)

# 讓 3D box 比例正確
ax.set_box_aspect((args.x_lim, args.y_lim, args.z_lim))

first_filename = os.path.basename(files[0])
title = ax.set_title(f"SPH 3D Sod Shock Tube - {first_filename}")

# 視角，可自行調整
ax.view_init(elev=20, azim=-60)

# 5. 更新動畫
def update(frame):
    filepath = files[frame]
    df = pd.read_csv(filepath)

    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    z = df["z"].to_numpy()
    rho = df["rho"].to_numpy()

    # 更新 3D scatter 位置
    sc._offsets3d = (x, y, z)

    # 更新顏色
    sc.set_array(rho)

    current_filename = os.path.basename(filepath)
    title.set_text(f"SPH 3D Sod Shock Tube - {current_filename}")

    return sc, title

# 6. 進度顯示
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 7. 建立動畫
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
    ani.save(output_file, writer="pillow", fps=args.fps,
             progress_callback=print_progress)

elif args.format == "mp4":
    ani.save(output_file, writer="ffmpeg", fps=args.fps,
             progress_callback=print_progress)

else:
    print("Wrong file format. Please use mp4 or gif.")
    exit(1)

print("\nDone!")

# python animate_3d_without_h.py -d ../sod3d_m0.001_x5.0_y1.0_z1.0_t5.0 -o sod3d_test -x 5.0 -y 1.0 -z 1.0 -f mp4 --fps 1