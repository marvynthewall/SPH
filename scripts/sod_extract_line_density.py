import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Initialize argument parser
parser = argparse.ArgumentParser(description="Plot SPH simulation density distribution")

# Add -f argument with a default value
parser.add_argument("-f", "--filename", type=str, default="output_0000.csv", 
                    help="CSV filename to read (default: output_0000.csv)")
parser.add_argument("-x", "--x", type=float, default=1.0, 
                    help="size of x direction (default: 1.0)")

# Parse arguments
args = parser.parse_args()

# Read data using the parsed filename
print(f"Reading file: {args.filename} ...")

def cubic_spline_kernel_2d(r, h):
    """
    2D Cubic Spline Kernel (支援單一數值或 numpy array 向量化運算)
    
    Parameters:
        r (float or np.ndarray): 粒子間距離
        h (float): Smoothing length
        
    Returns:
        tuple: (W, dWdr, dWdh)
    """
    # 確保 r 是 numpy array，方便進行向量化遮罩運算
    r = np.asarray(r, dtype=float)
    
    q = r / h
    norm = 40.0 / (7.0 * np.pi * h * h)
    norm_dW = norm / h
    
    # 初始化輸出陣列為 0 (這也同時處理了 q > 1.0 的情況)
    W = np.zeros_like(r)
    dWdr = np.zeros_like(r)
    
    # --- 區間 1: 0 <= q <= 0.5 ---
    mask1 = (q >= 0.0) & (q <= 0.5)
    q1 = q[mask1]
    W[mask1] = norm * (1.0 - 6.0 * q1**2 + 6.0 * q1**3)
    dWdr[mask1] = norm_dW * (-12.0 * q1 + 18.0 * q1**2)
    
    # --- 區間 2: 0.5 < q <= 1.0 ---
    mask2 = (q > 0.5) & (q <= 1.0)
    q2 = q[mask2]
    diff = 1.0 - q2
    W[mask2] = norm * 2.0 * diff**3
    dWdr[mask2] = norm_dW * (-6.0 * diff**2)
    
    # --- 計算 dWdh ---
    # 因為 W 和 dWdr 在 q > 1.0 時皆為 0，這個公式可以直接套用於整個陣列
    dWdh = (-2.0 / h) * W - q * dWdr
    
    return W


# 1. 讀取資料
df = pd.read_csv(args.filename)
x_part = df['x'].values
y_part = df['y'].values
m_part = df['m'].values
h = df['h'].values[0] # 假設 h 是常數

# 2. 定義抽樣線 (例如 y = 0.5, x 從 0 到 1 抽 500 個點)
y_target = 0.5
x_sample = np.linspace(0, args.x, 5000)
rho_sample = np.zeros_like(x_sample)

# 3. 計算密度 (向量化運算)
for i, x_s in enumerate(x_sample):
    # 計算該抽樣點到所有粒子的距離
    r = np.sqrt((x_part - x_s)**2 + (y_part - y_target)**2)

    # 篩選在 h 範圍內的粒子
    mask = r < h

    if np.any(mask):
        # 計算 Kernel 並加總質量
        w = cubic_spline_kernel_2d(r[mask], h)
        rho_sample[i] = np.sum(m_part[mask] * w)

# 4. 畫圖
plt.plot(x_sample, rho_sample, label=f'Density at y={y_target}')
plt.xlabel('x')
plt.ylabel('Density')
plt.show()
