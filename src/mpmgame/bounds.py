"""Chapter-2 bound helpers.

All bounds here are implemented as *numerical experiment-side checks*.
They should not be interpreted as theorem proofs.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .fmp import ContractSystem, compute_attack_map, compute_gamma, q_strength_metric


@dataclass
class BoundResult:
    value: float
    bound_type: str
    assumptions: str
    exact_for_static: bool


def bound_upper_full(system: ContractSystem, Pbar: np.ndarray, Cbar: np.ndarray | None = None) -> BoundResult:
    """Upper bound helper for full-threat comparisons.

    Implemented numerically as ||Gamma * Pbar||_2 for static systems (Q=0 setting).
    """
    Gamma = compute_gamma(system.G, system.alpha)
    n = system.n
    Q0 = np.zeros((n, n))
    M = compute_attack_map(Gamma, Q0, Pbar, Cbar=Cbar)
    value = float(np.linalg.svd(M, compute_uv=False)[0])
    return BoundResult(
        value=value,
        bound_type="upper_full",
        assumptions="Full attacker; static proxy in Q=0 realization; compatible Cbar/Pbar model required.",
        exact_for_static=True,
    )

# Gabriel Perry's lower bound
def bound_lower_single_link(system: ContractSystem, Q: np.ndarray, Pbar: np.ndarray) -> BoundResult:
    """Lower bound helper under Cbar=I single-link assumptions.

    V >= ||Pbar||_max / (n^(3/2) * (1 + max_i sigma_max(Q_i)) * max_i{1, ||G_i||_max}).
    Static simplification uses matrix max norms.
    """
    n = system.n
    pbar_max = float(np.max(np.abs(Pbar)))
    qterm = 1.0 + q_strength_metric(Q)
    gterm = max(1.0, float(np.max(np.abs(system.G))))
    denom = (n ** 1.5) * qterm * gterm
    value = pbar_max / denom if denom > 0 else 0.0
    return BoundResult(
        value=float(value),
        bound_type="lower_single_link",
        assumptions="Cbar=I; fixed single-link attack model; static simplification.",
        exact_for_static=False,
    )

# Gavin Glenn's lower bound
def bound_lower_full(system: ContractSystem, Pbar: np.ndarray) -> BoundResult:
    """Lower bound helper for full-threat case.

    V >= sigma_min(Pbar) / (1 + ||G alpha||).
    """
    svals = np.linalg.svd(np.asarray(Pbar, dtype=float), compute_uv=False)
    sigma_min = float(np.min(svals)) if svals.size else 0.0
    ga = np.asarray(system.G, dtype=float) @ np.asarray(system.alpha, dtype=float)
    ga_norm = float(np.linalg.svd(ga, compute_uv=False)[0])
    value = sigma_min / (1.0 + ga_norm)
    return BoundResult(
        value=float(value),
        bound_type="lower_full",
        assumptions="Full attacker; static simplification; compatible Pbar assumptions required.",
        exact_for_static=True,
    )


def check_bound_direction(measured: float, bound: BoundResult, direction: str, tol: float = 1e-9) -> bool:
    if direction == "upper":
        return measured <= bound.value + tol
    if direction == "lower":
        return measured + tol >= bound.value
    raise ValueError("direction must be 'upper' or 'lower'")
