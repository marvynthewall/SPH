import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# 1. 找出所有的輸出檔案並排序
file_pattern = "output_*.csv"
files = sorted(glob.glob(file_pattern))

if not files:
    print(f"Error: No files found matching '{file_pattern}'")
    exit(1)

print(f"Found {len(files)} files. Starting analysis...")

# 儲存計算結果的列表
times = []
Eky_list = []
Ek_list = []
U_list = []
Etotal_list = []
Px_list = []
Py_list = []

dt_output = 0.01  # 從 kh_2d.c 得知的輸出間隔

for idx, f in enumerate(files):
    df = pd.read_csv(f)
    
    # 時間 (根據輸出間隔推算)
    t = idx * dt_output
    times.append(t)
    
    # 質量、速度、內能
    m = df['m']
    vx = df['vx']
    vy = df['vy']
    u = df['u']
    
    # 計算 Y 方向動能
    Eky = np.sum(0.5 * m * vy**2)
    Eky_list.append(Eky)
    
    # 計算 X 方向動能與總動能
    Ekx = np.sum(0.5 * m * vx**2)
    Ek = Ekx + Eky
    Ek_list.append(Ek)
    
    # 計算總內能
    U = np.sum(m * u)
    U_list.append(U)
    
    # 計算系統總能量
    Etotal = Ek + U
    Etotal_list.append(Etotal)
    
    # 計算系統總動量
    Px = np.sum(m * vx)
    Py = np.sum(m * vy)
    Px_list.append(Px)
    Py_list.append(Py)
    
    if idx % 20 == 0:
        print(f"Processed {f} (t = {t:.2f})")

print("Analysis complete. Generating plots...")

# 2. 繪製圖表
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# --- Plot 1: KH 不穩定性成長率 (Y方向動能的對數圖) ---
axs[0].plot(times, np.log(Eky_list), 'b-', linewidth=2)
axs[0].set_title("Kelvin-Helmholtz Instability Growth Rate\n(Linear phase shows constant slope)", fontsize=14)
axs[0].set_xlabel("Time (t)")
axs[0].set_ylabel("ln(E_{k,y})")
axs[0].grid(True, linestyle='--', alpha=0.7)

# --- Plot 2: 系統能量守恆檢查 ---
# 為了看清楚變化，可以將能量減去初始值來觀察變化量，或者直接畫絕對值
axs[1].plot(times, Ek_list, 'g-', label='Kinetic Energy ($E_k$)')
axs[1].plot(times, U_list, 'r-', label='Internal Energy ($U$)')
axs[1].plot(times, Etotal_list, 'k-', linewidth=2, label='Total Energy ($E_{total}$)')
axs[1].set_title("Energy Conservation Check", fontsize=14)
axs[1].set_xlabel("Time (t)")
axs[1].set_ylabel("Energy")
axs[1].legend()
axs[1].grid(True, linestyle='--', alpha=0.7)

# --- Plot 3: 系統動量守恆檢查 ---
axs[2].plot(times, Px_list, 'm-', label='Total X-Momentum ($P_x$)')
axs[2].plot(times, Py_list, 'c-', label='Total Y-Momentum ($P_y$)')
axs[2].set_title("Momentum Conservation Check", fontsize=14)
axs[2].set_xlabel("Time (t)")
axs[2].set_ylabel("Momentum")
axs[2].legend()
axs[2].grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
output_image = "kh_analysis.png"
plt.savefig(output_image, dpi=150)
print(f"Plot saved to {output_image}")
