# ============================================================
# SPH Project Makefile
# ============================================================

# ----------------------------
# Phony targets
# ----------------------------
.PHONY: all clean run animate test

# ----------------------------
# Compiler settings
# ----------------------------
OMP ?= 0
GPU ?= 0

ifeq ($(GPU),1)
    # GPU compilation settings
	CC   = gcc
	NVCC = nvcc
	CFLAGS   = -O3 $(OPENFLAG) -Wall -Iinclude/
	NVCCFLAGS = -O3 -Iinclude/ -Xcompiler "-Wall"
	LDFLAGS  = -lm -L/usr/lib/x86_64-linux-gnu -lcudart 

else
    # CPU Openmp compilation settings
	ifeq ($(OMP),1)
		CC = gcc-14
		OPENFLAG = -fopenmp
	# CPU normal compilation settings
	else
		CC = gcc
		OPENFLAG =
	endif
	CFLAGS   = -O3 $(OPENFLAG) -Wall -Iinclude/
	LDFLAGS  = -lm $(OPENFLAG)
endif




# ----------------------------
# Linker flags
# ----------------------------
LDFLAGS = -lm $(OPENFLAG)


# ----------------------------
# Directories
# ----------------------------
BUILD_DIR = build
BIN_DIR   = bin


# ----------------------------
# Object files
# ----------------------------
OBJS = $(BUILD_DIR)/sph_system.o \
       $(BUILD_DIR)/kernel.o     \
       $(BUILD_DIR)/density.o    \
       $(BUILD_DIR)/force.o      \
       $(BUILD_DIR)/init.o       \
       $(BUILD_DIR)/integrator.o \
       $(BUILD_DIR)/io.o


# For GPU compilation
# ifeq ($(GPU),1)
#     OBJS += $(BUILD_DIR)/cuda_kernels.o
# endif


# ----------------------------
# Default target
# ----------------------------
all: $(BIN_DIR)/sod_2d
	@if [ "$(GPU)" = "1" ]; then \
		echo "[CUDA GPU Acceleration Enabled]"; \
	elif [ "$(OMP)" = "1" ]; then \
		echo "[OpenMP Parallel Acceleration Enabled]"; \
	else \
		echo "[Standard CPU Compilation]"; \
	fi

# ----------------------------
# Create directories
# ----------------------------
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)


# ----------------------------
# Main program
# ----------------------------
$(BIN_DIR)/sod_2d: $(OBJS) $(BUILD_DIR)/sod_2d.o | $(BIN_DIR)
	$(CC) ${CFLAGS} -o $@ $(OBJS) $(BUILD_DIR)/sod_2d.o $(LDFLAGS)


# ----------------------------
# Compile core source files
# ----------------------------
$(BUILD_DIR)/sph_system.o: src/sph_system.c include/sph_system.h | $(BUILD_DIR)
	$(CC) ${CFLAGS} -c src/sph_system.c -o $@

$(BUILD_DIR)/kernel.o: src/kernel.c include/kernel.h include/sph_system.h | $(BUILD_DIR)
	$(CC) ${CFLAGS} -c src/kernel.c -o $@

$(BUILD_DIR)/density.o: src/density.c include/density.h include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c src/density.c -o $@

$(BUILD_DIR)/force.o: src/force.c include/force.h include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c src/force.c -o $@

$(BUILD_DIR)/init.o: src/init.c include/init.h include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c src/init.c -o $@

$(BUILD_DIR)/integrator.o: src/integrator.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c src/integrator.c -o $@

$(BUILD_DIR)/io.o: src/io.c include/io.h include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c src/io.c -o $@


# For GPU compilation
$(BUILD_DIR)/cuda_kernels.o: src/cuda_kernels.cu include/sph_system.h | $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) -c src/cuda_kernels.cu -o $@


# ----------------------------
# Compile main example file
# ----------------------------
$(BUILD_DIR)/sod_2d.o: examples/sod_2d.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c examples/sod_2d.c -o $@

# ----------------------------
# Run main program
# ----------------------------
run: $(BIN_DIR)/sod_2d
	./$(BIN_DIR)/sod_2d
	cd $(BIN_DIR) && ./sod_2d

# ----------------------------
# Animation
# ----------------------------
animate: run
	python3 scripts/animate_2d.py

# ----------------------------
# Tests
# ----------------------------
test: tests/test_density tests/test_force tests/test_init tests/test_kernel tests/test_integrator

tests/test_density: $(BUILD_DIR)/test_density.o $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_force: $(BUILD_DIR)/test_force.o $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_init: $(BUILD_DIR)/test_init.o $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_kernel: $(BUILD_DIR)/test_kernel.o $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tests/test_integrator: $(BUILD_DIR)/test_integrator.o $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# ----------------------------
# Compile test files
# ----------------------------
$(BUILD_DIR)/test_density.o: tests/test_density.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c tests/test_density.c -o $@

$(BUILD_DIR)/test_force.o: tests/test_force.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c tests/test_force.c -o $@

$(BUILD_DIR)/test_init.o: tests/test_init.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c tests/test_init.c -o $@

$(BUILD_DIR)/test_kernel.o: tests/test_kernel.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c tests/test_kernel.c -o $@

$(BUILD_DIR)/test_integrator.o: tests/test_integrator.c include/sph_system.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c tests/test_integrator.c -o $@


# ----------------------------
# Clean
# ----------------------------
clean:
	rm -rf $(BUILD_DIR)      \
	       $(BIN_DIR)/sod_2d \
	       tests/test_kernel \
	       tests/test_density \
	       tests/test_force \
	       tests/test_init \
           tests/test_integrator
