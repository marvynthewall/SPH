import os
import argparse
import re
import matplotlib.pyplot as plt

def parse_time(filepath):
    """Extract Total Execution Time from a log file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        for line in f:
            if "Total Execution Time:" in line:
                match = re.search(r"Total Execution Time:\s*([0-9.]+)\s*seconds", line)
                if match:
                    return float(match.group(1))
    return None

def main():
    parser = argparse.ArgumentParser(description="Parse and plot SPH Performance Benchmarks")
    parser.add_argument("--logdir", type=str, default="benchmark_logs", help="Directory containing log txt files")
    parser.add_argument("-o", "--output", type=str, default="sod_performance.png", help="Output plot filename")
    args = parser.parse_args()

    cpu_time = parse_time(os.path.join(args.logdir, "cpu_log.txt"))
    gpu_time = parse_time(os.path.join(args.logdir, "gpu_log.txt"))

    threads_list = [1, 2, 4, 8, 16, 32]
    plot_threads = []
    omp_cell = []
    omp_nocell = []

    for th in threads_list:
        tc = parse_time(os.path.join(args.logdir, f"omp_{th}_log.txt"))
        tn = parse_time(os.path.join(args.logdir, f"omp_nocell_{th}_log.txt"))
        if tc is not None:
            plot_threads.append(th)
            omp_cell.append(tc)
            omp_nocell.append(tn)

    if not plot_threads and cpu_time is None and gpu_time is None:
        print(f"Error: No valid log files found in {args.logdir}")
        return

    # Base time for speedup calculation
    base_time = cpu_time if cpu_time is not None else (omp_cell[0] if omp_cell else 1.0)
    
    print("=== Performance Benchmark Results ===")
    if cpu_time: print(f"CPU (1 Core): {cpu_time:.3f}s")
    if gpu_time: print(f"GPU (CUDA): {gpu_time:.3f}s")
    for i, th in enumerate(plot_threads):
        tc = omp_cell[i]
        tn = omp_nocell[i]
        print(f"OMP ({th} Thr): Cell={tc:.3f}s" + (f", NO-CELL={tn:.3f}s" if tn else ""))

    # Plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Execution Time
    if omp_cell:
        ax1.plot(plot_threads, omp_cell, color='#1f77b4', marker='s', linestyle='-', linewidth=2, label="OMP (Cell-List, $\mathcal{O}(N)$)")
    
    valid_nocell_idx = [i for i, v in enumerate(omp_nocell) if v is not None]
    if valid_nocell_idx:
        x_nocell = [plot_threads[i] for i in valid_nocell_idx]
        y_nocell = [omp_nocell[i] for i in valid_nocell_idx]
        ax1.plot(x_nocell, y_nocell, color='#ff7f0e', marker='^', linestyle='-', linewidth=2, label="OMP (NO Cell-List, $\mathcal{O}(N^2)$)")

    if cpu_time:
        ax1.plot([1], [cpu_time], color='black', marker='*', markersize=14, label="CPU (1 Core Base)")
    if gpu_time:
        ax1.axhline(y=gpu_time, color='purple', linestyle=':', linewidth=2, label="GPU (CUDA)")

    ax1.set_xlabel('Number of Threads', fontsize=14)
    ax1.set_ylabel('Total Execution Time (seconds)', fontsize=14, color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.set_title('SPH Acceleration Performance: Cell-List vs Direct $O(N^2)$', fontsize=16)

    # Secondary Y-axis for Speedup
    ax2 = ax1.twinx()
    if omp_cell:
        speedup_cell = [base_time / t for t in omp_cell]
        ax2.plot(plot_threads, speedup_cell, color='#d62728', marker='o', linestyle='--', linewidth=2, label="Speedup (Cell-List)")
        for x, y in zip(plot_threads, speedup_cell):
            ax2.text(x, y - (max(speedup_cell)*0.05), f'{y:.1f}x', ha='center', va='top', color='#d62728', fontweight='bold', fontsize=10)

    if valid_nocell_idx:
        speedup_nocell = [base_time / y for y in y_nocell]
        ax2.plot(x_nocell, speedup_nocell, color='#8c564b', marker='v', linestyle='--', linewidth=2, label="Speedup (NO-CELL)")

    if cpu_time:
        ax2.plot([1], [1.0], color='black', marker='*')

    if gpu_time:
        speedup_gpu = base_time / gpu_time
        ax2.axhline(y=speedup_gpu, color='purple', linestyle=':', linewidth=2)
        ax2.text(plot_threads[-1] if plot_threads else 1, speedup_gpu + 0.1, f'GPU: {speedup_gpu:.1f}x', 
                 color='purple', fontweight='bold', ha='right', va='bottom')

    ax2.set_ylabel('Speedup (relative to CPU/Base)', fontsize=14, color='#d62728')
    ax2.tick_params(axis='y', labelcolor='#d62728')
    
    # Combine legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', bbox_to_anchor=(1.05, 1))

    if plot_threads:
        plt.xticks(plot_threads)
    
    fig.tight_layout()
    plt.savefig(args.output, dpi=200)
    print(f"Performance plot saved to {args.output}")

if __name__ == "__main__":
    main()
