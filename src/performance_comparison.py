"""Plot the CPU-time versus error comparison reported in the Master's thesis.

The numerical values are stored in results/comparison_results.csv. The thesis
figure is preserved at figures/cpu_vs_error.png and is not overwritten; this
script writes figures/cpu_vs_error_reproduced.png from the archived table.
This script does not invent or recompute new convergence studies; it reproduces
only the thesis-level CPU/error comparison.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "comparison_results.csv"
FIG = ROOT / "figures" / "cpu_vs_error_reproduced.png"


def main():
    df = pd.read_csv(DATA)
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.loglog(df["Euler_Error"], df["Euler_CPU_s"], "s-", label="Euler")
    ax.loglog(df["RK2_Error"], df["RK2_CPU_s"], "d-", label="RK2")
    ax.loglog(df["BDF_Error"], df["BDF_CPU_s"], "o-", label="BDF")

    ax.set_xlabel("Error (RMSE)")
    ax.set_ylabel("CPU Time (s)")
    ax.set_title("Comparison of Methods: CPU Time vs Error (log-log scale)")
    ax.grid(True, which="both", linestyle="--", alpha=0.45)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG, dpi=180)
    plt.close(fig)
    print(f"Wrote {FIG}")


if __name__ == "__main__":
    main()
