"""Finite zero-sum game solvers."""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog


def expected_utility(attacker_mix: np.ndarray, defender_mix: np.ndarray, payoff: np.ndarray) -> float:
    a = np.asarray(attacker_mix, dtype=float)
    d = np.asarray(defender_mix, dtype=float)
    A = np.asarray(payoff, dtype=float)
    return float(a @ A @ d)


def solve_zero_sum_game(payoff: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve a finite zero-sum matrix game.

    Attacker maximizes value, defender minimizes.
    """
    A = np.asarray(payoff, dtype=float)
    m, n = A.shape

    # Attacker LP: maximize v s.t. A^T p >= v, p>=0, sum p=1.
    # Convert to minimize -v with variables [p_1..p_m, v]
    c = np.zeros(m + 1)
    c[-1] = -1.0
    A_ub = np.zeros((n, m + 1))
    b_ub = np.zeros(n)
    for j in range(n):
        A_ub[j, :m] = -A[:, j]
        A_ub[j, -1] = 1.0
    A_eq = np.zeros((1, m + 1))
    A_eq[0, :m] = 1.0
    b_eq = np.array([1.0])
    bounds = [(0.0, None)] * m + [(None, None)]
    res_a = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if not res_a.success:
        raise RuntimeError(f"attacker LP failed: {res_a.message}")
    p = res_a.x[:m]
    v = res_a.x[-1]

    # Defender LP: minimize w s.t. A q <= w, q>=0, sum q=1.
    c2 = np.zeros(n + 1)
    c2[-1] = 1.0
    A_ub2 = np.zeros((m, n + 1))
    b_ub2 = np.zeros(m)
    for i in range(m):
        A_ub2[i, :n] = A[i, :]
        A_ub2[i, -1] = -1.0
    A_eq2 = np.zeros((1, n + 1))
    A_eq2[0, :n] = 1.0
    b_eq2 = np.array([1.0])
    bounds2 = [(0.0, None)] * n + [(None, None)]
    res_d = linprog(c2, A_ub=A_ub2, b_ub=b_ub2, A_eq=A_eq2, b_eq=b_eq2, bounds=bounds2, method="highs")
    if not res_d.success:
        raise RuntimeError(f"defender LP failed: {res_d.message}")
    q = res_d.x[:n]

    p /= p.sum()
    q /= q.sum()
    value = float(0.5 * (v + res_d.x[-1]))
    return p, q, value
