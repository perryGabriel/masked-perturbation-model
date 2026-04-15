"""LP-inspired static surrogate for chapter-2 design.

This is a finite-dimensional relaxation, not an exact solver for the
infinite-dimensional formulation in the notes.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import linprog

from .fmp import ContractSystem, access_matrix, compute_gamma, vulnerability_full, vulnerability_single_link, well_posed


@dataclass
class LPDesignResult:
    success: bool
    message: str
    Q: np.ndarray
    theta: np.ndarray
    surrogate_value: float
    true_vulnerability: float


def lp_relaxation_design(
    system: ContractSystem,
    Q_linearization: np.ndarray | None = None,
    access_model: str = "w2",
    threat_model: str = "single_link",
    g_hat: float = 0.2,
    structure_mask: np.ndarray | None = None,
    cond_max: float = 1e6,
) -> LPDesignResult:
    """Solve an LP on linearized static map M ≈ Gamma (I + Q) Pbar."""
    n = system.n
    Qk = np.zeros((n, n)) if Q_linearization is None else np.asarray(Q_linearization, dtype=float)

    Gamma = compute_gamma(system.G, system.alpha, cond_max=cond_max)
    Pbar = access_matrix(system, Qk, access_model)

    idx = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if structure_mask is not None and not bool(structure_mask[i, j]):
                continue
            idx.append((i, j))
    m = len(idx)

    M0 = Gamma @ Pbar
    Bmats = []
    for (i, j) in idx:
        E = np.zeros((n, n))
        E[i, j] = 1.0
        Bmats.append(Gamma @ (E @ Pbar))

    # Variable z = [theta_1..theta_m, t]
    c = np.zeros(m + 1)
    c[-1] = 1.0

    rows = []
    rhs = []
    # Entrywise sampled constraints for surrogate |M_ij| <= t
    # (single-link exact for static; full-case uses conservative entrywise surrogate)
    for a in range(n):
        for b in range(n):
            coeff = np.array([Bmats[k][a, b] for k in range(m)], dtype=float)
            base = float(M0[a, b])
            rows.append(np.r_[coeff, -1.0])
            rhs.append(-base)
            rows.append(np.r_[-coeff, -1.0])
            rhs.append(base)

    A_ub = np.array(rows, dtype=float) if rows else np.zeros((0, m + 1))
    b_ub = np.array(rhs, dtype=float) if rhs else np.zeros(0)

    bounds = [(-g_hat, g_hat) for _ in range(m)] + [(0.0, None)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res.success:
        return LPDesignResult(False, res.message, Qk.copy(), np.zeros(m), np.inf, np.inf)

    theta = res.x[:m]
    Q = np.zeros((n, n))
    for k, (i, j) in enumerate(idx):
        Q[i, j] = theta[k]

    if not well_posed(system, Q, cond_max=cond_max):
        return LPDesignResult(False, "Designed Q failed well-posedness check", Qk.copy(), theta, float(res.x[-1]), np.inf)

    if threat_model == "single_link":
        true_v = vulnerability_single_link(system, Q, access_matrix(system, Q, access_model)).value
    elif threat_model == "full":
        true_v = vulnerability_full(system, Q, access_matrix(system, Q, access_model)).value
    else:
        raise ValueError("threat_model must be 'single_link' or 'full'")

    return LPDesignResult(True, res.message, Q, theta, float(res.x[-1]), float(true_v))
