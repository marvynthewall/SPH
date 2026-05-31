# SPH

本專案是一個以 **C 語言**撰寫的 Smoothed Particle Hydrodynamics（SPH）模擬程式。  
專案架構拆分為多個模組，包含 SPH 系統資料結構、kernel function、密度計算、作用力計算、初始化、時間積分與輸出功能。

---

## 專案架構
```text
SPH/
├── Makefile
├── include/
├── src/
├── examples/
├── tests/
├── scripts/
└── bin/
```
---
## 編譯

### 純 CPU（預設）
```bash
make cpu
# 產生 bin/sod_2d_cpu
```

### OpenMP 多核心平行
```bash
make omp
# 產生 bin/sod_2d_omp
```

### CUDA GPU 加速
```bash
make gpu
# 產生 bin/sod_2d_gpu
```

### 三種全部編譯
```bash
make cpu omp gpu
```

---
## 執行

### 執行模擬
```bash
make run_cpu   # 純 CPU
make run_omp   # OpenMP
make run_gpu   # GPU
```

或直接執行：
```bash
./bin/sod_2d_cpu
./bin/sod_2d_omp
./bin/sod_2d_gpu
```

### 若需要測試：
全部測試：
```bash
make test
```
單一檔案測試，以 `test_density.c` 為例：
```bash
make tests/test_density
```


產生的測試執行檔會在 tests 內，目前測試程式包含：

```text
tests/test_density
tests/test_force
tests/test_init
tests/test_kernel
```
對應的測試原始碼位於：
```text
tests/test_density.c
tests/test_force.c
tests/test_init.c
tests/test_kernel.c
```


---
### 產生動畫
```bash
make animate
```


### 清除編譯檔案
```bash
make clean
```

會刪除所有執行檔，以及編譯相關檔案：
```text
build
bin/sod_2d
bin/sod_2d_cpu
bin/sod_2d_omp
bin/sod_2d_gpu
bin/output*.csv
tests/test_density
tests/test_force
tests/test_init
tests/test_kernel
tests/test_integrator
```