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

parser.add_argument("-x", "--x_lim", type=float, default=15.0,
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

# 4. 建立 3D figure - 將畫布比例調得更扁（12x3），逼迫長條形震波管展開
fig = plt.figure(figsize=(20, 10))

# 讓 3D 繪圖軸極致外擴，上下左右甚至故意超出畫布（溢出），直接剪裁掉 3D 引擎自帶的死白邊
ax = fig.add_axes([-0.18, -0.25, 1.36, 1.45], projection="3d")

sc = ax.scatter(
    df["x"],
    df["y"],
    df["z"],
    c=df["rho"],
    s=args.point_size,
    cmap="viridis",
    vmin=0.0,
    vmax=1.0,
    depthshade=True
)

# 🌟 終極特技 1：改用「正交投影」
# 消除透視造成的遠小近大，讓 15 單位長的 X 軸以完全均勻的比例橫跨螢幕
ax.set_proj_type('ortho')

# 🌟 終極特技 2：強力鏡頭拉近（放大三倍）
# Matplotlib 3.6+ 支援 set_zoom()。設為 2.8 ~ 3.0 可以直接把本體放大三倍！
# 為了相容舊版本，若不支援則自動 fallback 到修改相機距離 ax.dist（數字越小鏡頭越近，預設為 10）
try:
    ax.set_zoom(2.9)
except AttributeError:
    ax.dist = 3.5

# 🌟 Colorbar 絕對定位：水平放置於左下角，高度與邊距依新的畫布比例微調
cax = fig.add_axes([0.05, 0.08, 0.20, 0.05]) 
cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
cbar.set_label("Density", labelpad=5)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.set_xlim(0, args.x_lim)
ax.set_ylim(0, args.y_lim)
ax.set_zlim(0, args.z_lim)

# 讓 3D box 比例正確 (15:1:1)
ax.set_box_aspect((args.x_lim, args.y_lim, args.z_lim))

# 🌟 全域標題置頂防切：y 軸鎖定在 0.95 確保大圖時字體清晰完整
first_filename = os.path.basename(files[0])
title = fig.suptitle(f"SPH 3D Sod Shock Tube - {first_filename}", y=0.95, fontsize=12)

# 視角調整：維持你習慣的角度
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
