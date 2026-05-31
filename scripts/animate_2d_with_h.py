import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import glob
import argparse
import os

# 1. 初始化參數解析
parser = argparse.ArgumentParser(description="Create SPH animation with particles and smoothing length circles.")
parser.add_argument("-d", "--dir", type=str, required=True, 
                    help="Target directory containing output_*.csv files")
parser.add_argument("-o", "--outputfilename", type=str, default="sod_animation_with_h", 
                    help="Output gif filename (default: sod_animation_with_h)")
parser.add_argument("-x", "--x_lim", type=float, default=5.0, 
                    help="Size of x direction limit (default: 5.0)")
parser.add_argument("-f", "--format", type=str, default="gif", 
                    help="File format, mp4 or gif (default: gif)")
args = parser.parse_args()

# 2. 尋找並排序 CSV 檔案
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}!")
    exit(1)

print(f"Found {len(files)} files. Initializing plot...")

# 3. 建立畫布與初始設定
fig, ax = plt.subplots(figsize=(15, 3))

# 讀取第一幀資料
df = pd.read_csv(files[0])

# --- 繪製初始的 h 圓形 (背景層) ---
circles = [Circle((x, y), r) for x, y, r in zip(df['x'], df['y'], df['h'])]
patch_col = PatchCollection(circles, alpha=0.15, facecolor='gray', edgecolor='none', zorder=1)
ax.add_collection(patch_col)

# --- 繪製初始的粒子中心 (前景層) ---
# s=5 是點的大小，可以根據你的視覺需求調整
sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=5, zorder=2, vmin=0.0, vmax=1.2, cmap='viridis')

# 圖表格式設定
fig.colorbar(sc, ax=ax, label="Density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(0, args.x_lim)
ax.set_ylim(0, 1.0)

# 強制 x, y 軸等比例，避免圓形變成橢圓
ax.set_aspect('equal')

first_filename = os.path.basename(files[0])
title = ax.set_title(f"SPH 2D Sod Shock Tube - {first_filename}")

# 4. 定義動畫更新函數
def update(frame):
    global patch_col # 宣告 global 以便在此函數內覆寫全域的 patch_col
    
    filepath = files[frame]
    df = pd.read_csv(filepath)
    
    # (1) 更新散佈圖 (粒子中心位置與密度顏色)
    sc.set_offsets(df[["x", "y"]].values)
    sc.set_array(df["rho"].values)
    
    # (2) 更新 h 圓形集合
    # 由於 PatchCollection 不易直接修改每個圓的屬性，最穩健的做法是移除舊的並加入新的
    patch_col.remove() 
    new_circles = [Circle((x, y), r) for x, y, r in zip(df['x'], df['y'], df['h'])]
    patch_col = PatchCollection(new_circles, alpha=0.15, facecolor='gray', edgecolor='none', zorder=1)
    ax.add_collection(patch_col)
    
    # (3) 更新標題
    current_filename = os.path.basename(filepath)
    title.set_text(f"SPH 2D Sod Shock Tube - {current_filename}")
    
    return sc, patch_col, title

# 5. 進度條回呼函數
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rSaving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 6. 產生與儲存動畫
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

output_file = args.outputfilename + "." + args.format
print(f"Saving animation to {output_file} ...")
if args.format == 'gif':
    ani.save(output_file, writer='pillow', fps=10, progress_callback=print_progress)
elif args.format == 'mp4':
    ani.save(output_file, writer='ffmpeg', fps=10, progress_callback=print_progress)
else:
    print("wrong file format")
print("\nDone!")
