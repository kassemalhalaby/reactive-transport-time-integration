"""Explicit RK2 (midpoint) for the finite-volume method-of-lines system.

This is the scheme reported in the thesis:

    k1 = f(C^n)
    k2 = f(C^n + (dt/2) k1)
    C^{n+1} = C^n + dt k2

with the same fixed step dt = 4 s used for implicit Euler, so that the two
fixed-step integrators are compared at identical temporal resolution.

Being explicit, the step is bounded by a stability condition. For this operator
the diffusive bound dominates: dt <= min R(C) dx^2 / (2 D), and since R(C) >= 1
the conservative limit is dx^2 / (2 D). At ne = 500 that is about 173 s, so
dt = 4 s is well inside the stable range; at ne = 4000 it falls to about 2.7 s
and the step becomes stability-bound rather than accuracy-bound.
"""
import numpy as np

L = 1.0
T = 0.4 * 24 * 3600
V = 1.0 / (24 * 3600)
D = 1e-3 / (24 * 3600)
eps = 0.4
rho_b = 1.59
Kf = 0.126
nf = 0.7
dt = 4.0
CLEFT = 0.0
CRIGHT = 0.0


def R(C):
    return 1.0 + (rho_b / eps) * Kf * nf * np.maximum(C, 1e-12) ** (nf - 1.0)


def initial_condition(x):
    C = np.full_like(x, 1e-4, dtype=float)
    m = (0.15 <= x) & (x < 0.2); C[m] = (x[m] - 0.15) / 0.05
    m = (0.2 <= x) & (x < 0.3); C[m] = 1.0
    m = (0.3 <= x) & (x < 0.35); C[m] = 1.0 - (x[m] - 0.3) / 0.05
    return C


def rhs(t, C, dx):
    dCdt = np.zeros_like(C)
    dCdt[0] = (-V * (C[1] - CLEFT) / (2 * dx)
               + D * (C[1] - 2 * C[0] + CLEFT) / dx**2) / R(C[0])
    dCdt[1:-1] = (-V * (C[2:] - C[:-2]) / (2 * dx)
                  + D * (C[2:] - 2 * C[1:-1] + C[:-2]) / dx**2) / R(C[1:-1])
    dCdt[-1] = (-V * (CRIGHT - C[-2]) / (2 * dx)
                + D * (CRIGHT - 2 * C[-1] + C[-2]) / dx**2) / R(C[-1])
    return dCdt


def stability_limit(ne):
    """Conservative explicit step bound dx^2 / (2 D), in seconds."""
    return (L / ne) ** 2 / (2.0 * D)


def solve(ne=500, step=dt):
    dx = L / ne
    xc = np.linspace(dx / 2.0, L - dx / 2.0, ne)
    C = initial_condition(xc)

    t = 0.0
    for _ in range(int(round(T / step))):
        k1 = rhs(t, C, dx)
        k2 = rhs(t + 0.5 * step, C + 0.5 * step * k1, dx)
        C = C + step * k2
        t += step
    return xc, C


if __name__ == "__main__":
    ne = 500
    print(f"ne = {ne}, dt = {dt} s, stability limit = {stability_limit(ne):.1f} s")
    x, C = solve(ne)
    np.savetxt("rk2_result.txt", np.column_stack((x, C)), header="x C")
