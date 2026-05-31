# ============================================================
# SPH Project Makefile
# ============================================================

.PHONY: all cpu omp gpu clean run run_cpu run_omp run_gpu animate test

# ----------------------------
# Directories
# ----------------------------
BUILD_CPU = build/cpu
BUILD_OMP = build/omp
BUILD_GPU = build/gpu
BIN_DIR   = bin

# ----------------------------
# Compiler settings
# ----------------------------
CC_CPU    = gcc
CC_OMP    = gcc-14
NVCC      = nvcc

CFLAGS_CPU  = -O3 -Wall -Iinclude/
CFLAGS_OMP  = -O3 -Wall -Iinclude/ -fopenmp
NVCCFLAGS   = -O3 -Iinclude/ -Xcompiler "-Wall"

LDFLAGS_CPU = -lm
LDFLAGS_OMP = -lm -fopenmp
LDFLAGS_GPU = -lm -L/usr/lib/x86_64-linux-gnu -lcudart

# ----------------------------
# Object files
# ----------------------------
SRCS_COMMON = sph_system kernel init io

OBJS_CPU = $(foreach s, $(SRCS_COMMON), $(BUILD_CPU)/$(s).o) \
           $(BUILD_CPU)/density.o \
           $(BUILD_CPU)/force.o   \
           $(BUILD_CPU)/integrator.o

OBJS_OMP = $(foreach s, $(SRCS_COMMON), $(BUILD_OMP)/$(s).o) \
           $(BUILD_OMP)/density.o \
           $(BUILD_OMP)/force.o   \
           $(BUILD_OMP)/integrator.o

OBJS_GPU = $(foreach s, $(SRCS_COMMON), $(BUILD_GPU)/$(s).o) \
           $(BUILD_GPU)/density.o \
           $(BUILD_GPU)/force.o   \
           $(BUILD_GPU)/integrator.o

# ----------------------------
# Default target
# ----------------------------
all: cpu kh_2d

cpu: $(BIN_DIR)/sod_2d_cpu
	@echo "[Standard CPU Compilation]"

omp: $(BIN_DIR)/sod_2d_omp
	@echo "[OpenMP Parallel Acceleration Enabled]"

gpu: $(BIN_DIR)/sod_2d_gpu
	@echo "[CUDA GPU Acceleration Enabled]"

# ----------------------------
# Create directories
# ----------------------------
$(BUILD_CPU):
	mkdir -p $(BUILD_CPU)

$(BUILD_OMP):
	mkdir -p $(BUILD_OMP)

$(BUILD_GPU):
	mkdir -p $(BUILD_GPU)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# ----------------------------
# Link
# ----------------------------
$(BIN_DIR)/sod_2d_cpu: $(OBJS_CPU) $(BUILD_CPU)/sod_2d.o | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

$(BIN_DIR)/sod_2d_omp: $(OBJS_OMP) $(BUILD_OMP)/sod_2d.o | $(BIN_DIR)
	$(CC_OMP) $(CFLAGS_OMP) -o $@ $^ $(LDFLAGS_OMP)

$(BIN_DIR)/sod_2d_gpu: $(OBJS_GPU) $(BUILD_GPU)/sod_2d.o | $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(LDFLAGS_GPU)

kh_2d: $(BIN_DIR)/kh_2d

$(BIN_DIR)/kh_2d: $(OBJS_CPU) $(BUILD_CPU)/kh_2d.o | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

# ----------------------------
# Compile examples
# ----------------------------
$(BUILD_CPU)/sod_2d.o: examples/sod_2d.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/kh_2d.o: examples/kh_2d.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_OMP)/sod_2d.o: examples/sod_2d.c | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_GPU)/sod_2d.o: examples/sod_2d.c | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ----------------------------
# Compile common sources (CPU)
# ----------------------------
$(BUILD_CPU)/sph_system.o: src/sph_system.c include/sph_system.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/kernel.o: src/kernel.c include/kernel.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/init.o: src/init.c include/init.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/io.o: src/io.c include/io.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/density.o: src/density.c include/density.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/force.o: src/force.c include/force.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/integrator.o: src/integrator.c include/integrator.h | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

# ----------------------------
# Compile common sources (OMP)
# ----------------------------
$(BUILD_OMP)/sph_system.o: src/sph_system.c include/sph_system.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/kernel.o: src/kernel.c include/kernel.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/init.o: src/init.c include/init.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/io.o: src/io.c include/io.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/density.o: src/density.c include/density.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/force.o: src/force.c include/force.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_OMP)/integrator.o: src/integrator.c include/integrator.h | $(BUILD_OMP)
	$(CC_OMP) $(CFLAGS_OMP) -c $< -o $@

$(BUILD_DIR)/kh_2d.o: examples/kh_2d.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c examples/kh_2d.c -o $@

# ----------------------------
# Compile common sources (GPU)
# ----------------------------
$(BUILD_GPU)/sph_system.o: src/sph_system.c include/sph_system.h | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/kernel.o: src/kernel.c include/kernel.h | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/init.o: src/init.c include/init.h | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/io.o: src/io.c include/io.h | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/density.o: src/density.cu include/density.cuh | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/force.o: src/force.cu include/force.cuh | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_GPU)/integrator.o: src/integrator.cu include/integrator.cuh | $(BUILD_GPU)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ----------------------------
# Run
# ----------------------------
run_cpu: cpu
	./$(BIN_DIR)/sod_2d_cpu

run_omp: omp
	./$(BIN_DIR)/sod_2d_omp

run_gpu: gpu
	./$(BIN_DIR)/sod_2d_gpu

run: run_cpu

# ----------------------------
# Animation
# ----------------------------
animate: run
	python3 scripts/animate_2d.py

# ----------------------------
# Tests
# ----------------------------
test: tests/test_density tests/test_force tests/test_init tests/test_kernel tests/test_integrator

tests/test_density: $(BUILD_CPU)/test_density.o $(OBJS_CPU) | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

tests/test_force: $(BUILD_CPU)/test_force.o $(OBJS_CPU) | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

tests/test_init: $(BUILD_CPU)/test_init.o $(OBJS_CPU) | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

tests/test_kernel: $(BUILD_CPU)/test_kernel.o $(OBJS_CPU) | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

tests/test_integrator: $(BUILD_CPU)/test_integrator.o $(OBJS_CPU) | $(BIN_DIR)
	$(CC_CPU) $(CFLAGS_CPU) -o $@ $^ $(LDFLAGS_CPU)

$(BUILD_CPU)/test_density.o: tests/test_density.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/test_force.o: tests/test_force.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/test_init.o: tests/test_init.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/test_kernel.o: tests/test_kernel.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

$(BUILD_CPU)/test_integrator.o: tests/test_integrator.c | $(BUILD_CPU)
	$(CC_CPU) $(CFLAGS_CPU) -c $< -o $@

# ----------------------------
# Clean
# ----------------------------
clean:
	rm -rf build/          \
	       bin/sod_2d      \
	       bin/sod_2d_cpu  \
	       bin/sod_2d_omp  \
	       bin/sod_2d_gpu  \
	       bin/kh_2d       \
	       bin/output*.csv \
	       tests/test_density    \
	       tests/test_force      \
	       tests/test_init       \
	       tests/test_kernel     \
	       tests/test_integrator
