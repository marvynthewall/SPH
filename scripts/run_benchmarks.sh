#!/bin/bash
# ==============================================================================
# Script to run SPH Sod 2D benchmarks across CPU, OMP, and GPU
# It saves the output logs to a text file for parsing by python.
# ==============================================================================

# Benchmark parameters
MASS="0.001"
X="15.0"
TEND="4.0"

# Create log directory
mkdir -p benchmark_logs

echo "======================================"
echo " Starting Performance Benchmarks"
echo "======================================"

# 1. Single Core CPU
echo "Running Single Core CPU..."
if [ -f "bin/sod_2d_cpu" ]; then
    ./bin/sod_2d_cpu -m $MASS -x $X -t $TEND -o benchmark_logs/cpu_out > benchmark_logs/cpu_log.txt
    echo "  -> Done"
else
    echo "  -> Skipping (bin/sod_2d_cpu not found)"
fi

# 2. OpenMP (Various Threads) - Cell-List
echo "Running OpenMP (Cell-List)..."
if [ -f "bin/sod_2d_omp" ]; then
    for threads in 1 2 4 8 16; do
        echo "  - $threads threads (Cell-List)"
        ./bin/sod_2d_omp -m $MASS -x $X -t $TEND -th $threads -o benchmark_logs/omp_${threads}_out > benchmark_logs/omp_${threads}_log.txt
    done
    
    echo "Running OpenMP (NO Cell-List)..."
    for threads in 1 2 4 8 16; do
        echo "  - $threads threads (NO Cell-List)"
        ./bin/sod_2d_omp -m $MASS -x $X -t $TEND -th $threads -nocell -o benchmark_logs/omp_nocell_${threads}_out > benchmark_logs/omp_nocell_${threads}_log.txt
    done
    echo "  -> Done"
else
    echo "  -> Skipping (bin/sod_2d_omp not found)"
fi

# 3. GPU
echo "Running GPU..."
if [ -f "bin/sod_2d_gpu" ]; then
    ./bin/sod_2d_gpu -m $MASS -x $X -t $TEND -o benchmark_logs/gpu_out > benchmark_logs/gpu_log.txt
    echo "  -> Done"
else
    echo "  -> Skipping (bin/sod_2d_gpu not found)"
fi

echo "======================================"
echo " Benchmarks Complete!"
echo " Please run: python scripts/benchmark_performance.py --logdir benchmark_logs"
echo "======================================"
