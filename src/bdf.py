"""Variable-order BDF experiment used in the Master's thesis.

The spatial operator is the same cell-centred finite-volume central stencil used
for the other time integrators.  The ODE system is passed to SciPy's variable-
order BDF integrator.  An analytical Jacobian is supplied, following the
nonlinear derivatives used in the thesis implementation.
"""
import numpy as np
from scipy.integrate import solve_ivp

L = 1.0
T = 0.4 * 24 * 3600
V = 1.0 / (24 * 3600)
D = 1e-3 / (24 * 3600)
eps = 0.4
rho_b = 1.59
Kf = 0.126
nf = 0.7
CLEFT = 0.0
CRIGHT = 0.0


def R(C):
    C = np.maximum(C, 1e-12)
    return 1.0 + (rho_b / eps) * Kf * nf * C ** (nf - 1.0)


def dR_dC(C):
    C = np.maximum(C, 1e-12)
    return (rho_b / eps) * Kf * nf * (nf - 1.0) * C ** (nf - 2.0)


def initial_condition(x):
    C = np.full_like(x, 1e-4, dtype=float)
    m = (0.15 <= x) & (x < 0.2)
    C[m] = (x[m] - 0.15) / 0.05
    m = (0.2 <= x) & (x < 0.3)
    C[m] = 1.0
    m = (0.3 <= x) & (x < 0.35)
    C[m] = 1.0 - (x[m] - 0.3) / 0.05
    return C


def rhs(t, C, dx):
    out = np.zeros_like(C)
    out[0] = (-V * (C[1] - CLEFT) / (2 * dx)
              + D * (C[1] - 2 * C[0] + CLEFT) / dx**2) / R(C[0])
    out[1:-1] = (-V * (C[2:] - C[:-2]) / (2 * dx)
                 + D * (C[2:] - 2 * C[1:-1] + C[:-2]) / dx**2) / R(C[1:-1])
    out[-1] = (-V * (CRIGHT - C[-2]) / (2 * dx)
               + D * (CRIGHT - 2 * C[-1] + C[-2]) / dx**2) / R(C[-1])
    return out


def jacobian(t, C, dx):
    """Analytical Jacobian of the method-of-lines right-hand side.

    Tridiagonal: each cell couples only to its two neighbours. The diagonal
    carries the quotient-rule term -R'(C) q(C) / R(C)^2, where q is the spatial
    operator; the boundary rows use the same central advective closure as rhs().
    """
    n = len(C)
    J = np.zeros((n, n))

    # Left boundary row
    r0 = R(C[0])
    dr0 = dR_dC(C[0])
    q0 = (-V * (C[1] - CLEFT) / (2 * dx)
          + D * (C[1] - 2 * C[0] + CLEFT) / dx**2)
    J[0, 0] = (-dr0 * q0 / r0**2) - 2 * D / (r0 * dx**2)
    J[0, 1] = (-V / (2 * dx) + D / dx**2) / r0

    # Interior rows
    for i in range(1, n - 1):
        ri = R(C[i])
        dri = dR_dC(C[i])
        qi = (-V * (C[i + 1] - C[i - 1]) / (2 * dx)
              + D * (C[i + 1] - 2 * C[i] + C[i - 1]) / dx**2)
        J[i, i - 1] = (V / (2 * dx) + D / dx**2) / ri
        J[i, i] = -dri * qi / ri**2 - 2 * D / (ri * dx**2)
        J[i, i + 1] = (-V / (2 * dx) + D / dx**2) / ri

    # Right boundary row
    rN = R(C[-1])
    drN = dR_dC(C[-1])
    qN = (-V * (CRIGHT - C[-2]) / (2 * dx)
          + D * (CRIGHT - 2 * C[-1] + C[-2]) / dx**2)
    J[-1, -2] = (V / (2 * dx) + D / dx**2) / rN
    J[-1, -1] = -drN * qN / rN**2 - 2 * D / (rN * dx**2)

    return J


def solve(ne=500, rtol=1e-6, atol=1e-9):
    dx = L / ne
    xc = np.linspace(dx / 2, L - dx / 2, ne)
    C0 = initial_condition(xc)

    sol = solve_ivp(
        lambda t, y: rhs(t, y, dx),
        (0.0, T),
        C0,
        method="BDF",
        jac=lambda t, y: jacobian(t, y, dx),
        rtol=rtol,
        atol=atol,
    )
    return xc, sol.y[:, -1], sol


if __name__ == "__main__":
    x, C, sol = solve()
    np.savetxt("bdf_result.txt", np.column_stack((x, C)), header="x C")
    print(f"BDF accepted steps: {len(sol.t) - 1}")
