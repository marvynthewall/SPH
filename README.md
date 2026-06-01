# SPH 模擬程式

以 C 語言撰寫的 Smoothed Particle Hydrodynamics（SPH）模擬程式，支援純 CPU、OpenMP 多核心平行以及 CUDA GPU 加速。
目前包含三種主要模擬案例：

- 2D Sod shock tube
- 3D Sod shock tube
- 2D Kelvin-Helmholtz instability

---

## 專案架構

```text
SPH/
├── Makefile
├── include/        # 標頭檔
├── src/            # 核心原始碼
├── examples/       # 範例模擬
│   ├── sod_2d.c
│   ├── sod_3d.c
│   └── kh_2d.c
├── tests/          # 單元測試
├── scripts/        # Python 後處理腳本
└── bin/            # 編譯後執行檔
```

---

## 環境需求

### 基本需求

- GCC（建議 ≥ 9）
- Make
- Python 3（選用，用於後處理與動畫）

### OpenMP 版本

- GCC-14 或支援 OpenMP 的 GCC 編譯器

### GPU 版本

- CUDA Toolkit
- NVIDIA GPU
- `nvcc` 編譯器

### Python 後處理需求

```bash
pip install numpy pandas matplotlib
```

若要輸出 mp4 動畫，需安裝 `ffmpeg`。

---

# Makefile 使用方式

本專案使用 `Makefile` 管理編譯與執行流程。

主要編譯規則如下：

```text
make             # 編譯全部模擬與全部後端
make sod_2d      # 編譯 2D Sod shock tube
make sod_3d      # 編譯 3D Sod shock tube
make kh_2d       # 編譯 2D Kelvin-Helmholtz

make cpu         # 編譯所有 CPU 版本
make omp         # 編譯所有 OpenMP 版本
make gpu         # 編譯所有 GPU 版本
```

主要執行規則如下：

```text
make run_sod2d_cpu
make run_sod2d_omp
make run_sod2d_gpu

make run_sod3d_cpu
make run_sod3d_omp
make run_sod3d_gpu

make run_kh_cpu
make run_kh_omp
make run_kh_gpu
```

---

# 編譯

## 1. 編譯全部模擬與全部後端

```bash
make
```

等同於：

```bash
make sod_2d sod_3d kh_2d
```

會產生：

```text
bin/sod_2d_cpu
bin/sod_2d_omp
bin/sod_2d_gpu

bin/sod_3d_cpu
bin/sod_3d_omp
bin/sod_3d_gpu

bin/kh_2d_cpu
bin/kh_2d_omp
bin/kh_2d_gpu
```

---

## 2. 只編譯 2D Sod shock tube

```bash
make sod_2d
```

會產生：

```text
bin/sod_2d_cpu
bin/sod_2d_omp
bin/sod_2d_gpu
```

---

## 3. 只編譯 3D Sod shock tube

```bash
make sod_3d
```

會產生：

```text
bin/sod_3d_cpu
bin/sod_3d_omp
bin/sod_3d_gpu
```

---

## 4. 只編譯 2D Kelvin-Helmholtz instability

```bash
make kh_2d
```

會產生：

```text
bin/kh_2d_cpu
bin/kh_2d_omp
bin/kh_2d_gpu
```

---

## 5. 依照後端編譯

### 編譯所有 CPU 版本

```bash
make cpu
```

會產生：

```text
bin/sod_2d_cpu
bin/sod_3d_cpu
bin/kh_2d_cpu
```

---

### 編譯所有 OpenMP 版本

```bash
make omp
```

會產生：

```text
bin/sod_2d_omp
bin/sod_3d_omp
bin/kh_2d_omp
```

---

### 編譯所有 GPU 版本

```bash
make gpu
```

會產生：

```text
bin/sod_2d_gpu
bin/sod_3d_gpu
bin/kh_2d_gpu
```

---

# 執行

## 1. 執行 2D Sod shock tube

### CPU 版本

```bash
make run_sod2d_cpu
```

或直接執行：

```bash
./bin/sod_2d_cpu
```

---

### OpenMP 版本

```bash
make run_sod2d_omp
```

或直接執行：

```bash
./bin/sod_2d_omp
```

---

### GPU 版本

```bash
make run_sod2d_gpu
```

或直接執行：

```bash
./bin/sod_2d_gpu
```

---

## 2. 執行 3D Sod shock tube

### CPU 版本

```bash
make run_sod3d_cpu
```

或直接執行：

```bash
./bin/sod_3d_cpu
```

---

### OpenMP 版本

```bash
make run_sod3d_omp
```

或直接執行：

```bash
./bin/sod_3d_omp
```

---

### GPU 版本

```bash
make run_sod3d_gpu
```

或直接執行：

```bash
./bin/sod_3d_gpu
```

---

## 3. 執行 2D Kelvin-Helmholtz instability

### CPU 版本

```bash
make run_kh_cpu
```

或直接執行：

```bash
./bin/kh_2d_cpu
```

---

### OpenMP 版本

```bash
make run_kh_omp
```

或直接執行：

```bash
./bin/kh_2d_omp
```

---

### GPU 版本

```bash
make run_kh_gpu
```

或直接執行：

```bash
./bin/kh_2d_gpu
```

---

## 4. 預設執行

```bash
make run
```

等同於：

```bash
make run_sod2d_cpu
```

---

# 常用工作流程

## 只跑 2D Sod CPU

```bash
make clean
make run_sod2d_cpu
```

---

## 只跑 2D Sod OpenMP

```bash
make clean
make run_sod2d_omp
```

---

## 只跑 3D Sod CPU

```bash
make clean
make run_sod3d_cpu
```

---

## 只跑 3D Sod OpenMP

```bash
make clean
make run_sod3d_omp
```

---

## 只跑 Kelvin-Helmholtz CPU

```bash
make clean
make run_kh_cpu
```

---

## 只編譯所有 CPU 版本，不執行

```bash
make clean
make cpu
```

會編譯：

```text
bin/sod_2d_cpu
bin/sod_3d_cpu
bin/kh_2d_cpu
```

---

## 只編譯所有 OpenMP 版本，不執行

```bash
make clean
make omp
```

會編譯：

```text
bin/sod_2d_omp
bin/sod_3d_omp
bin/kh_2d_omp
```

---

## 只編譯所有 GPU 版本，不執行

```bash
make clean
make gpu
```

會編譯：

```text
bin/sod_2d_gpu
bin/sod_3d_gpu
bin/kh_2d_gpu
```

---

# 輸出資料

模擬執行後會產生一個輸出資料夾，內含：

```text
output_0000.csv
output_0001.csv
output_0002.csv
...
time_log.csv
```

每個 `output_XXXX.csv` 儲存該時間點的粒子資訊，例如：

```text
id,x,y,z,vx,vy,vz,ax,ay,az,m,rho,P,u,h,cs
```

其中：

| 欄位 | 意義 |
|---|---|
| `x, y, z` | 粒子位置 |
| `vx, vy, vz` | 粒子速度 |
| `ax, ay, az` | 粒子加速度 |
| `m` | 粒子質量 |
| `rho` | 密度 |
| `P` | 壓力 |
| `u` | specific internal energy |
| `h` | smoothing length |
| `cs` | sound speed |

---

# 動畫與後處理

## 1. 使用 Makefile 產生 2D Sod 動畫

```bash
make animate
```

此指令會先執行：

```bash
make run_sod2d_cpu
```

接著執行：

```bash
python3 scripts/animate_2d_with_h.py -d sod_m0.001_x5.0_t5.0
```

---

## 2. 手動產生 2D 動畫

```bash
python3 scripts/animate_2d_with_h.py \
    -d sod_m0.001_x5.0_t5.0 \
    -o sod_2d_animation \
    -x 5.0 \
    -f mp4
```

若要輸出 GIF：

```bash
python3 scripts/animate_2d_with_h.py \
    -d sod_m0.001_x5.0_t5.0 \
    -o sod_2d_animation \
    -x 5.0 \
    -f gif
```

---

## 3. 手動產生 3D 動畫

若 3D Sod 輸出資料夾為：

```text
sod3d_m0.001_x5.0_y1.0_z1.0_t5.0
```

可以執行：

```bash
python3 scripts/animate_3d_with_h.py \
    -d sod3d_m0.001_x5.0_y1.0_z1.0_t5.0 \
    -o sod3d_animation \
    -x 5.0 \
    -y 1.0 \
    -z 1.0 \
    -f mp4
```

如果目前人在 `scripts/` 資料夾內，且輸出資料夾位於上一層，則使用：

```bash
python3 animate_3d_with_h.py \
    -d ../sod3d_m0.001_x5.0_y1.0_z1.0_t5.0 \
    -o sod3d_animation \
    -x 5.0 \
    -y 1.0 \
    -z 1.0 \
    -f mp4
```

---

# 單元測試

## 執行全部測試

```bash
make test
```

---

## 執行單一測試

### Density 測試

```bash
make tests/test_density
./tests/test_density
```

### Force 測試

```bash
make tests/test_force
./tests/test_force
```

### Init 測試

```bash
make tests/test_init
./tests/test_init
```

### Kernel 測試

```bash
make tests/test_kernel
./tests/test_kernel
```

### Integrator 測試

```bash
make tests/test_integrator
./tests/test_integrator
```

---

# 清除編譯檔案

```bash
make clean
```

會刪除：

```text
build/

bin/sod_2d_cpu
bin/sod_2d_omp
bin/sod_2d_gpu

bin/sod_3d_cpu
bin/sod_3d_omp
bin/sod_3d_gpu

bin/kh_2d_cpu
bin/kh_2d_omp
bin/kh_2d_gpu

bin/output*.csv

tests/test_density
tests/test_force
tests/test_init
tests/test_kernel
tests/test_integrator
```

---

# 指令總表

| 指令 | 功能 |
|---|---|
| `make` | 編譯全部模擬與全部後端 |
| `make cpu` | 編譯所有 CPU 版本 |
| `make omp` | 編譯所有 OpenMP 版本 |
| `make gpu` | 編譯所有 GPU 版本 |
| `make sod_2d` | 編譯 2D Sod 的 CPU / OMP / GPU 版本 |
| `make sod_3d` | 編譯 3D Sod 的 CPU / OMP / GPU 版本 |
| `make kh_2d` | 編譯 KH 2D 的 CPU / OMP / GPU 版本 |
| `make run` | 預設執行 2D Sod CPU 版本 |
| `make run_sod2d_cpu` | 執行 2D Sod CPU 版本 |
| `make run_sod2d_omp` | 執行 2D Sod OpenMP 版本 |
| `make run_sod2d_gpu` | 執行 2D Sod GPU 版本 |
| `make run_sod3d_cpu` | 執行 3D Sod CPU 版本 |
| `make run_sod3d_omp` | 執行 3D Sod OpenMP 版本 |
| `make run_sod3d_gpu` | 執行 3D Sod GPU 版本 |
| `make run_kh_cpu` | 執行 KH 2D CPU 版本 |
| `make run_kh_omp` | 執行 KH 2D OpenMP 版本 |
| `make run_kh_gpu` | 執行 KH 2D GPU 版本 |
| `make animate` | 執行 2D Sod CPU 並產生動畫 |
| `make test` | 編譯並執行全部單元測試 |
| `make clean` | 清除編譯產物與執行檔 |
