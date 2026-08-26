"""Recreate the unified terminal-profile comparison.

Runs the three integrators on the same mesh and superposes the terminal
profiles with the initial condition and the digitized new-ELLAM reference.

The thesis figure is preserved untouched at figures/comparison_overlay.png.
This script writes figures/comparison_overlay_reproduced.png so that the
archived figure and a freshly computed one can be compared side by side.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import implicit_euler_central as euler
import rk2
import bdf

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "new_ellam_reference.csv"
OUT = ROOT / "figures" / "comparison_overlay_reproduced.png"

NE = 500


def main():
    dx = 1.0 / NE
    xc = np.linspace(dx / 2.0, 1.0 - dx / 2.0, NE)
    C0 = euler.initial_condition(xc)

    print(f"ne = {NE}, dt = {euler.dt} s")
    print("  implicit Euler ...")
    _, C_euler = euler.solve(ne=NE)
    print("  explicit RK2 ...")
    _, C_rk2 = rk2.solve(ne=NE)
    print("  variable-order BDF ...")
    _, C_bdf, sol = bdf.solve(ne=NE)
    print(f"  BDF accepted steps: {len(sol.t) - 1}")

    print("\nmaximum pointwise differences at T:")
    print(f"  |Euler - BDF| = {np.max(np.abs(C_euler - C_bdf)):.2e}")
    print(f"  |RK2   - BDF| = {np.max(np.abs(C_rk2 - C_bdf)):.2e}")

    ref = np.loadtxt(REF, delimiter=",", skiprows=1)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(xc, C0, "k--", lw=2, label="Initial Condition")
    ax.plot(xc, C_euler, "-", color="green", lw=2, label="Euler")
    ax.plot(xc, C_bdf, "-.", color="red", lw=2, label="BDF")
    ax.plot(xc, C_rk2, "--", color="blue", lw=2, label="RK2")
    ax.plot(ref[:, 0], ref[:, 1], "--", color="tab:cyan", lw=1.6,
            label="new-ELLAM curve")

    ax.set_xlabel("Distance x (m)")
    ax.set_ylabel("Concentration C")
    ax.set_title("Comparison of Solutions")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT, dpi=180)
    plt.close(fig)
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
