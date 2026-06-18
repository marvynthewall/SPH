#!/bin/bash
# ==============================================================================
# Script to run convergence tests (varying resolutions) for Sod Shock Tube
# Compares Grid 1D, SPH 1D, SPH 2D, SPH 3D.
# ==============================================================================

# Simulation physical time
TEND="1.0"
X="15.0"

# Mass values for SPH (lower mass = higher resolution)
MASSES=("0.01" "0.001") #這裡直接改mass

echo "======================================"
echo " Starting Convergence Benchmarks"
echo "======================================"

mkdir -p convergence_data

for m in "${MASSES[@]}"; do
    # Dynamically calculate Grid nx = X / mass (aligned with SPH 1D particle number)
    n=$(awk -v x="$X" -v m="$m" 'BEGIN {print int(x/m)}')
    
    echo "--- Running (mass=$m, Grid nx=$n) ---"
    
    # 1. Grid 1D
    if [ -f "grid_solver/bin/grid_sod" ]; then
        out_dir="convergence_data/grid_n${n}"
        echo "Running Grid 1D (nx=$n)..."
        ./grid_solver/bin/grid_sod -n $n -o $out_dir
    fi
    
    # 2. SPH 1D (using OMP)
    if [ -f "bin/sod_1d_omp" ]; then
        out_dir="convergence_data/sph1d_m${m}"
        echo "Running SPH 1D (m=$m)..."
        ./bin/sod_1d_omp -m $m -x $X -t $TEND -o $out_dir
    fi
    
    # 3. SPH 2D (using OMP)
    if [ -f "bin/sod_2d_omp" ]; then
        out_dir="convergence_data/sph2d_m${m}"
        echo "Running SPH 2D (m=$m)..."
        ./bin/sod_2d_omp -m $m -x $X -t $TEND -o $out_dir
    fi
    
    # 4. SPH 3D (using OMP)
    if [ -f "bin/sod_3d_omp" ]; then
        out_dir="convergence_data/sph3d_m${m}"
        echo "Running SPH 3D (m=$m)..."
        # 3D is extremely slow, maybe limit the mass for 3D or skip the highest resolution if needed.
        # But we will try running it.
        ./bin/sod_3d_omp -m $m -x $X -t $TEND -o $out_dir
    fi
done

echo "======================================"
echo " Convergence Benchmarks Complete!"
echo " Please run: python scripts/benchmark_convergence.py"
echo "======================================"
