# SPH_tommy

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

## OpenMP 設定

目前 OpenMP 預設是關閉的：

```makefile
# OPENFLAG=-fopenmp
```

如果需要啟用 OpenMP，可以將上面這行改成：

```makefile
OPENFLAG=-fopenmp
```

然後重新編譯：

```bash
make clean
make
```

---

## 基本工作流程
Makefile 會自動建立 `build/`，存放編譯所需檔案，不會 push 到 Github：

### 編譯 bin/sod_2d
```bash
make
```

### 編譯 bin/sod_2d，並且執行模擬
```bash
make run
```
等同於 make 之後，執行
```bash
./bin/sod_2d
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

### 清除編譯檔案
```bash
make clean
```

會刪除所有執行檔，以及編譯相關檔案：
```text
build/
bin/sod_2d
tests/test_kernel
tests/test_density
tests/test_force
tests/test_init
```