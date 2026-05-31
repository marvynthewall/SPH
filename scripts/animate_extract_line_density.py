import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os

# 1. 初始化參數解析
parser = argparse.ArgumentParser(description="Create 1D Line Density SPH animation for multiple Y slices.")
parser.add_argument("-d", "--dir", type=str, required=True, 
                    help="Target directory containing output_*.csv files")
parser.add_argument("-o", "--outputfilename", type=str, default="sod_line_density", 
                    help="Output filename without extension (default: sod_line_density)")
parser.add_argument("-x", "--x_lim", type=float, default=5.0, 
                    help="Size of x direction limit (default: 5.0)")
# 修改這裡：使用 nargs='+' 來接收多個數值，並給予預設值 [0.5]
parser.add_argument("-y", "--y_targets", type=float, nargs='+', default=[0.5], 
                    help="One or more Y coordinates for the sampling lines (e.g., -y 0.2 0.5 0.8)")
parser.add_argument("-f", "--format", type=str, default="gif", 
                    help="File format, mp4 or gif (default: gif)")
args = parser.parse_args()

# 2. Kernel 函數定義
def cubic_spline_kernel_2d(r, h):
    r = np.asarray(r, dtype=float)
    q = r / h
    norm = 40.0 / (7.0 * np.pi * h * h)
    
    W = np.zeros_like(r)
    
    # 區間 1: 0 <= q <= 0.5
    mask1 = (q >= 0.0) & (q <= 0.5)
    q1 = q[mask1]
    W[mask1] = norm * (1.0 - 6.0 * q1**2 + 6.0 * q1**3)
    
    # 區間 2: 0.5 < q <= 1.0
    mask2 = (q > 0.5) & (q <= 1.0)
    q2 = q[mask2]
    diff = 1.0 - q2
    W[mask2] = norm * 2.0 * diff**3
    
    return W

# 3. 尋找並排序 CSV 檔案
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}!")
    exit(1)

print(f"Found {len(files)} files. Target Y slices: {args.y_targets}")
print("Initializing plot...")

# 定義 1000 個 X 軸的抽樣點
x_sample = np.linspace(0, args.x_lim, 1000)

# 4. 建立畫布與初始設定
fig, ax = plt.subplots(figsize=(10, 6))

# 動態建立多條空折線
lines = []
for y_val in args.y_targets:
    # 每次呼叫 plot，Matplotlib 會自動分配不同顏色
    line, = ax.plot([], [], lw=2, alpha=0.5, label=f'y = {y_val}')
    lines.append(line)

# 圖表格式設定
ax.set_xlabel("x")
ax.set_ylabel("Density")
ax.set_xlim(0, args.x_lim)
ax.set_ylim(0, 1.2) 
ax.legend(loc="upper right")
ax.grid(True, linestyle='--', alpha=0.6)

first_filename = os.path.basename(files[0])
title = ax.set_title(f"1D Density Profile - {first_filename}")

# 5. 定義動畫更新函數
def update(frame):
    filepath = files[frame]
    df = pd.read_csv(filepath)
    
    x_part = df['x'].values
    y_part = df['y'].values
    m_part = df['m'].values
    h = df['h'].values[0] # 假設每一幀的 h 都是常數
    
    # 針對每一個指定的 y_target 進行計算與更新
    for idx, y_target in enumerate(args.y_targets):
        rho_sample = np.zeros_like(x_sample)
        
        for i, x_s in enumerate(x_sample):
            r = np.sqrt((x_part - x_s)**2 + (y_part - y_target)**2)
            mask = r < h
            if np.any(mask):
                w = cubic_spline_kernel_2d(r[mask], h)
                rho_sample[i] = np.sum(m_part[mask] * w)
                
        # 更新該條線的資料
        lines[idx].set_data(x_sample, rho_sample)
    
    current_filename = os.path.basename(filepath)
    title.set_text(f"1D Density Profile - {current_filename}")
    
    # 回傳所有更新的物件 (展開 list 並加上 title)
    return lines + [title]

# 6. 進度條回呼函數
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rComputing and saving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 7. 產生與儲存動畫
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

# 根據 format 參數決定輸出格式
if args.format == 'gif':
    output_file = args.outputfilename + ".gif"
    print(f"Saving animation to {output_file} ... (This will take a while!)")
    ani.save(output_file, writer='pillow', fps=10, progress_callback=print_progress)
elif args.format == 'mp4':
    output_file = args.outputfilename + ".mp4"
    print(f"Saving animation to {output_file} ... (This will take a while!)")
    ani.save(output_file, writer='ffmpeg', fps=10, progress_callback=print_progress)
else:
    print(f"Error: Unsupported format '{args.format}'. Please use 'gif' or 'mp4'.")

print("\nDone!")
