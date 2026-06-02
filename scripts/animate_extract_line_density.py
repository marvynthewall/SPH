import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import argparse
import os
import random

# 引入解析解函數
from sod_exact import get_sod_exact_rho

# 1. 初始化參數解析
parser = argparse.ArgumentParser(description="Create 1D Line Density SPH animation for multiple Y slices.")
parser.add_argument("-d", "--dir", type=str, required=True, 
                    help="Target directory containing output_*.csv and time_log.csv")
parser.add_argument("-o", "--outputfilename", type=str, default="sod_line_density", 
                    help="Output filename without extension (default: sod_line_density)")
parser.add_argument("-x", "--x_lim", type=float, default=5.0, 
                    help="Size of x direction limit (default: 5.0)")
parser.add_argument("-y", "--y_targets", type=float, nargs='+', 
                    help="One or more Y coordinates for the sampling lines (e.g., -y 0.2 0.5 0.8)")
parser.add_argument("-n", "--n_targets", type=int, default =5, 
                    help="Y coordinates random sample N lines (default: 5)")
parser.add_argument("-f", "--format", type=str, default="mp4", 
                    help="File format, mp4 or gif (default: mp4)")
parser.add_argument("-a", "--avg", type=bool, default=True, 
                    help="averaging the sample Y coordinates (default: True)")
args = parser.parse_args()

# 2. Kernel 函數定義
def cubic_spline_kernel_2d(r, h):
    r = np.asarray(r, dtype=float)
    q = r / h
    norm = 40.0 / (7.0 * np.pi * h * h)
    
    W = np.zeros_like(r)
    
    mask1 = (q >= 0.0) & (q <= 0.5)
    q1 = q[mask1]
    W[mask1] = norm * (1.0 - 6.0 * q1**2 + 6.0 * q1**3)
    
    mask2 = (q > 0.5) & (q <= 1.0)
    q2 = q[mask2]
    diff = 1.0 - q2
    W[mask2] = norm * 2.0 * diff**3
    
    return W

# 3. 尋找並排序 CSV 檔案，以及讀取時間紀錄
search_pattern = os.path.join(args.dir, "output_*.csv")
files = sorted(glob.glob(search_pattern))

if not files:
    print(f"No output files found in directory: {args.dir}!")
    exit(1)

# 讀取 time_log.csv
time_log_path = os.path.join(args.dir, "time_log.csv")
if not os.path.exists(time_log_path):
    print(f"Error: {time_log_path} not found in directory: {args.dir}!")
    exit(1)
time_log_df = pd.read_csv(time_log_path)

print(f"Found {len(files)} files. Target Y slices: {args.y_targets}")
print("Initializing plot...")

x_sample = np.linspace(0, args.x_lim, 5000)

# 4. 建立畫布與初始設定
fig, ax = plt.subplots(figsize=(10, 6))

if args.y_targets is not None:
    y_targets = args.y_targets
else:
    y_targets = [random.random() for _ in range(args.n_targets)]

print("Y samples at: ")
print(y_targets)


lines = []
if args.avg:
    line, = ax.plot([], [], lw=2, alpha=0.8, color='tab:blue', label='SPH Average Density')
    lines.append(line)
else:
    for y_val in y_targets:
        line, = ax.plot([], [], lw=2, alpha=0.5, label=f'SPH y = {y_val}')
        lines.append(line)

# 新增：建立一條解析解的線 (黑色虛線)
exact_line, = ax.plot([], [], 'k--', lw=1.5, label='Analytical Solution')

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
    current_filename = os.path.basename(filepath)
    
    # 從檔名解析出 step 數字 (例如 output_0010.csv -> 10)
    step_num = int(current_filename.replace("output_", "").replace(".csv", ""))
    
    # 從 time_log_df 找出對應的時間 t
    # 使用 .loc 確保能精準對應 frame 欄位，避免檔案順序錯亂導致時間對不起來
    current_t = time_log_df.loc[time_log_df['frame'] == step_num, 't'].values[0]
    
    x_part = df['x'].values
    y_part = df['y'].values
    m_part = df['m'].values
    h = df['h'].values[0] 
    
# 💡 根據 args.avg 走不同的密度計算與線條更新邏輯V
    if args.avg:
        # 建立一個累加密度的陣列
        rho_sum = np.zeros_like(x_sample)
        
        for y_target in y_targets:
            rho_sample = np.zeros_like(x_sample)
            for i, x_s in enumerate(x_sample):
                r = np.sqrt((x_part - x_s)**2 + (y_part - y_target)**2)
                mask = r < h
                if np.any(mask):
                    w = cubic_spline_kernel_2d(r[mask], h)
                    rho_sample[i] = np.sum(m_part[mask] * w)
            rho_sum += rho_sample
        
        # 取所有 y_target 的平均值，並更新到唯一的 lines[0]
        rho_avg = rho_sum / len(y_targets)
        lines[0].set_data(x_sample, rho_avg)
        
    else:
        # 原本的邏輯：每條 y_target 各自畫一條線
        for idx, y_target in enumerate(y_targets):
            rho_sample = np.zeros_like(x_sample)

            for i, x_s in enumerate(x_sample):
                r = np.sqrt((x_part - x_s)**2 + (y_part - y_target)**2)
                mask = r < h
                if np.any(mask):
                    w = cubic_spline_kernel_2d(r[mask], h)
                    rho_sample[i] = np.sum(m_part[mask] * w)

            lines[idx].set_data(x_sample, rho_sample)
        
    # 計算並更新解析解
    # x0 設定為管子的正中央 (args.x_lim / 2.0)
    rho_exact = get_sod_exact_rho(x_sample, current_t, x0=(args.x_lim / 2.0))
    exact_line.set_data(x_sample, rho_exact)
    
    # 更新標題，加入時間 t 顯示
    title.set_text(f"1D Density Profile - {current_filename} (t={current_t:.4f})")
    
    return lines + [exact_line, title]

# 6. 進度條回呼函數
def print_progress(current_frame, total_frames):
    percent = (current_frame + 1) / total_frames * 100
    print(f"\rComputing and saving frame {current_frame + 1}/{total_frames} ({percent:.1f}%)", end="")

# 7. 產生與儲存動畫
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)

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
