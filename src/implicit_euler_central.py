"""Implicit Euler with central finite-volume advection.

Core method used in the Master's thesis. The fully implicit nonlinear step is
solved with Newton-Raphson using the analytical tridiagonal Jacobian.
"""
import numpy as np

# Benchmark parameters from the thesis
L = 1.0
T = 0.4 * 24 * 3600
V = 1.0 / (24 * 3600)
D = 1e-3 / (24 * 3600)
eps = 0.4
rho_b = 1.59
Kf = 0.126
nf = 0.7
dt = 4.0
TOL = 1e-6
MAX_ITER = 20
CLEFT = 0.0
CRIGHT = 0.0


def R(C):
    return 1.0 + (rho_b / eps) * Kf * nf * np.maximum(C, 1e-12) ** (nf - 1.0)


def dR_dC(C):
    C = np.maximum(C, 1e-12)
    return (rho_b / eps) * Kf * nf * (nf - 1.0) * C ** (nf - 2.0)


def initial_condition(x):
    C = np.full_like(x, 1e-4, dtype=float)
    C[(0.15 <= x) & (x < 0.2)] = (x[(0.15 <= x) & (x < 0.2)] - 0.15) / 0.05
    C[(0.2 <= x) & (x < 0.3)] = 1.0
    C[(0.3 <= x) & (x < 0.35)] = 1.0 - (x[(0.3 <= x) & (x < 0.35)] - 0.3) / 0.05
    return C


def solve(ne=500):
    """One implicit Euler run.

    The linear system is tridiagonal, but it is solved here with a dense
    numpy.linalg.solve, which is the implementation used to produce the CPU
    times reported in the thesis and archived in results/comparison_results.csv.
    Keeping it dense is what makes those timings reproducible.
    """
    dx = L / ne
    xc = np.linspace(dx / 2.0, L - dx / 2.0, ne)
    Cold = initial_condition(xc)
    C = Cold.copy()
    ndt = int(T / dt)

    for _ in range(ndt):
        for _ in range(MAX_ITER):
            J = np.zeros((ne, ne))
            F = np.zeros(ne)

            # Left boundary row
            J[0, 0] = (dx / dt) * (R(C[0]) + (C[0] - Cold[0]) * dR_dC(C[0])) + 2 * D / dx
            J[0, 1] = V / 2 - D / dx
            F[0] = (dx / dt) * R(C[0]) * (C[0] - Cold[0]) \
                   + V * (C[1] - CLEFT) / 2 \
                   - (D / dx) * (CLEFT - 2 * C[0] + C[1])

            for i in range(1, ne - 1):
                J[i, i - 1] = -V / 2 - D / dx
                J[i, i] = (dx / dt) * (R(C[i]) + (C[i] - Cold[i]) * dR_dC(C[i])) + 2 * D / dx
                J[i, i + 1] = V / 2 - D / dx
                F[i] = (dx / dt) * R(C[i]) * (C[i] - Cold[i]) \
                       + V * (C[i + 1] - C[i - 1]) / 2 \
                       - (D / dx) * (C[i + 1] - 2 * C[i] + C[i - 1])

            # Right boundary row
            J[-1, -2] = -V / 2 - D / dx
            J[-1, -1] = (dx / dt) * (R(C[-1]) + (C[-1] - Cold[-1]) * dR_dC(C[-1])) + 2 * D / dx
            F[-1] = (dx / dt) * R(C[-1]) * (C[-1] - Cold[-1]) \
                    + V * (CRIGHT - C[-2]) / 2 \
                    - (D / dx) * (CRIGHT - 2 * C[-1] + C[-2])

            delta = np.linalg.solve(J, -F)
            C += delta
            if np.linalg.norm(delta, 2) < TOL:
                break
        Cold = C.copy()

    return xc, C


if __name__ == "__main__":
    x, C = solve()
    np.savetxt("implicit_euler_central_result.txt", np.column_stack((x, C)), header="x C")
