import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os

# 1. 初始化參數解析
parser = argparse.ArgumentParser(description="Create SPH animation with particles.")
parser.add_argument("-d", "--dir", type=str, required=True, 
                    help="Target directory containing output_*.csv files")
parser.add_argument("-o", "--outputfilename", type=str, default="sod_animation", 
                    help="Output filename without extension (default: sod_animation)")
parser.add_argument("-x", "--x_lim", type=float, default=15.0, 
                    help="Size of x direction limit (default: 15.0)")
parser.add_argument("-f", "--format", type=str, default="mp4", 
                    help="File format, mp4 or gif (default: mp4)")
args = parser.parse_args()

# 2. 尋找並排序 CSV 檔案
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}!")
    exit(1)

print(f"Found {len(files)} files. Initializing plot...")

# 3. 建立畫布與初始設定
# 調整畫布比例（例如 16:4），並加上 layout='constrained' 自動優化排版，減少周圍留白
fig, ax = plt.subplots(figsize=(16, 4), layout='constrained')

# 讀取第一幀資料
df = pd.read_csv(files[0])

N = len(df)
dynamic_s = np.interp(N, [500, 8000, 80000], [5.0, 2.0, 1.0])
# --- 繪製初始的粒子中心 (前景層) ---
# 【調整】將 s 從 5 縮小到 1 或 2，讓點點看起來更細緻
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=1.5, zorder=2, vmin=0.0, vmax=1.2, cmap='viridis')

# 圖表格式設定
# 【調整】使用 pad 控制 colorbar 與主圖的距離，shrink 控制它的長度（不讓它撐開上下留白）
cbar = fig.colorbar(sc, ax=ax, label="Density", pad=0.02, shrink=0.7)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(0, args.x_lim)
ax.set_ylim(0, 1.0)

# 強制 x, y 軸等比例，避免圓形變成橢圓
ax.set_aspect('equal')

# 【調整】將標題稍微往上抬，或利用 pad 增加間距，避免壓縮到主圖
first_filename = os.path.basename(files[0])
title = ax.set_title(f"SPH 2D Sod Shock Tube - {first_filename}", pad=10)

# 4. 定義動畫更新函數
def update(frame):
    filepath = files[frame]
    df = pd.read_csv(filepath)
    
    # 更新散佈圖 (粒子中心位置與密度顏色)
    sc.set_offsets(df[["x", "y"]].values)
    sc.set_array(df["rho"].values)
    
    # 更新標題
    current_filename = os.path.basename(filepath)
    title.set_text(f"SPH 2D Sod Shock Tube - {current_filename}")
    
    return sc, title

# 5. 進度條回呼函數
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 6. 產生與儲存動畫
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

# 根據 format 參數決定輸出格式
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
