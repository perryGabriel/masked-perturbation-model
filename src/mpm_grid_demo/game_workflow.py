"""Full ANDES/MPM finite-action game workflow helpers.

These utilities preserve the numerical workflow from the standalone
``mpm-andes-demo`` notebook while moving the reusable pieces into the original
``masked-perturbation-model`` repository.

The model-map object used for frequency-domain scans is the transfer-function
matrix ``M`` returned by :func:`mpmgame.lft.build_model_map`.  Static
attack-success checks intentionally use the equivalent closed-loop state matrix

    A_cl(Delta, mask) = A + B_w (Delta ∘ mask) C_r,

which matches the original ANDES notebook and avoids fragile MIMO
transfer-function conversion requirements in ``python-control``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import control
import numpy as np
import pandas as pd

import mpmgame as mpm


@dataclass(frozen=True)
class HinfApproxResult:
    norm: float
    critical_frequency: float
    frequency_grid: np.ndarray
    singular_values: np.ndarray
    critical_left_vector: np.ndarray | None = None
    critical_right_vector: np.ndarray | None = None


@dataclass
class StaticReducedGame:
    payoff: np.ndarray
    attacks: list[mpm.AttackAction]
    defenses: list[mpm.DefenseAction]


@dataclass
class StaticSuccessSets:
    attack_success: dict[str, set[str]]
    defense_success: dict[str, set[str]]


def frequency_grid(w_min: float, w_max: float, n_points: int) -> np.ndarray:
    """Return a log-spaced positive frequency grid."""
    w_min = float(w_min)
    w_max = float(w_max)
    n_points = int(n_points)
    if w_min <= 0 or w_max <= 0:
        raise ValueError("frequency_min and frequency_max must be positive")
    if w_max <= w_min:
        raise ValueError("frequency_max must be larger than frequency_min")
    if n_points < 2:
        raise ValueError("frequency_points must be at least 2")
    return np.logspace(np.log10(w_min), np.log10(w_max), n_points)


def freqresp_matrix(M: np.ndarray, omega: float) -> np.ndarray:
    """Evaluate a transfer-function matrix at ``s = j omega``."""
    M = np.asarray(M, dtype=object)
    if M.ndim != 2:
        raise ValueError("M must be a 2D transfer-function matrix")
    out = np.empty(M.shape, dtype=complex)
    s = 1j * float(omega)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            out[i, j] = complex(control.evalfr(M[i, j], s))
    return out


def approximate_hinf_norm(M: np.ndarray, omega_grid: Sequence[float]) -> HinfApproxResult:
    """Approximate ``||M||_∞`` by gridding the maximum singular value."""
    omega = np.asarray(omega_grid, dtype=float)
    sigmas = np.zeros_like(omega, dtype=float)
    best_idx = 0
    best_u: np.ndarray | None = None
    best_v: np.ndarray | None = None

    for k, w in enumerate(omega):
        Gjw = freqresp_matrix(M, float(w))
        U, S, Vh = np.linalg.svd(Gjw, full_matrices=False)
        sigmas[k] = float(S[0]) if S.size else 0.0
        if sigmas[k] >= sigmas[best_idx]:
            best_idx = k
            best_u = U[:, 0] if S.size else None
            best_v = Vh.conj().T[:, 0] if S.size else None

    return HinfApproxResult(
        norm=float(sigmas[best_idx]),
        critical_frequency=float(omega[best_idx]),
        frequency_grid=omega,
        singular_values=sigmas,
        critical_left_vector=best_u,
        critical_right_vector=best_v,
    )


def single_link_vulnerability_heatmap(M: np.ndarray, omega_grid: Sequence[float]) -> np.ndarray:
    """Return ``V[i,j] = max_ω |M_{j,i}(jω)|``.

    Rows index write/input channels and columns index read/output channels,
    matching the shape of the attack matrix ``Delta`` and defense mask ``∇``.
    """
    M = np.asarray(M, dtype=object)
    if M.ndim != 2:
        raise ValueError("M must be a 2D transfer-function matrix")
    n_reads, n_writes = M.shape
    V = np.zeros((n_writes, n_reads), dtype=float)
    for w in omega_grid:
        Gjw = freqresp_matrix(M, float(w))
        V = np.maximum(V, np.abs(Gjw).T)
    return V


def top_k_single_link_attacks(V: np.ndarray, k: int) -> list[tuple[int, int, float]]:
    """Return the top entries of ``V`` as ``(write_idx, read_idx, value)``."""
    V = np.asarray(V, dtype=float)
    if V.ndim != 2:
        raise ValueError("V must be a 2D array")
    entries = [
        (write_idx, read_idx, float(V[write_idx, read_idx]))
        for write_idx in range(V.shape[0])
        for read_idx in range(V.shape[1])
    ]
    entries.sort(key=lambda item: item[2], reverse=True)
    return entries[: int(k)]


def rank1_static_delta_from_hinf(
    hinf_result: HinfApproxResult,
    safety_factor: float = 1.05,
) -> np.ndarray:
    """Construct a real static rank-one ``Delta`` from peak singular vectors."""
    u = hinf_result.critical_left_vector
    v = hinf_result.critical_right_vector
    if u is None or v is None or hinf_result.norm <= 0:
        raise ValueError("cannot construct rank-one attack from a zero/empty Hinf result")
    return np.asarray(
        (float(safety_factor) / hinf_result.norm) * np.real(np.outer(v, np.conjugate(u))),
        dtype=float,
    )


def single_link_static_delta(
    write_index: int,
    read_index: int,
    vulnerability_value: float,
    safety_factor: float,
    n_writes: int,
    n_reads: int,
) -> np.ndarray:
    """Construct a sparse static ``Delta`` with one read→write link."""
    Delta = np.zeros((int(n_writes), int(n_reads)), dtype=float)
    gain = 0.0 if vulnerability_value <= 0 else float(safety_factor) / float(vulnerability_value)
    Delta[int(write_index), int(read_index)] = gain
    return Delta


def masked_static_delta(delta: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Apply an MPM defense mask to a static ``Delta``."""
    delta = np.asarray(delta, dtype=float)
    mask = np.asarray(mask, dtype=float)
    if delta.shape != mask.shape:
        raise ValueError(f"delta shape {delta.shape} and mask shape {mask.shape} must match")
    return delta * mask


def closed_loop_matrix_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    Delta: np.ndarray,
) -> np.ndarray:
    """Return ``A_cl = A + B_w Delta C_r`` for static output feedback."""
    A = np.asarray(A, dtype=float)
    B_w = np.asarray(B_w, dtype=float)
    C_r = np.asarray(C_r, dtype=float)
    Delta = np.asarray(Delta, dtype=float)
    return A + B_w @ Delta @ C_r


def closed_loop_eigs_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    Delta: np.ndarray,
) -> np.ndarray:
    """Return closed-loop eigenvalues for static output feedback."""
    return np.linalg.eigvals(closed_loop_matrix_static(A, B_w, C_r, Delta))


def is_destabilizing_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    delta: np.ndarray,
    mask: np.ndarray | None = None,
    pole_tol: float = 1e-8,
) -> bool:
    """True if ``A + B_w (delta ∘ mask) C_r`` is not strictly stable."""
    effective_delta = np.asarray(delta, dtype=float)
    if mask is not None:
        effective_delta = masked_static_delta(effective_delta, mask)
    eigs = closed_loop_eigs_static(A, B_w, C_r, effective_delta)
    return bool(np.max(np.real(eigs)) >= -float(pole_tol))


def compute_success_sets_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    attacks: Sequence[mpm.AttackAction],
    defenses: Sequence[mpm.DefenseAction],
    pole_tol: float = 1e-8,
) -> StaticSuccessSets:
    """Compute attack/defense success sets using static closed-loop eigenvalues."""
    attack_success: dict[str, set[str]] = {}
    defense_success: dict[str, set[str]] = {}

    for attack in attacks:
        attack_success[attack.label] = set()
        for defense in defenses:
            if is_destabilizing_static(A, B_w, C_r, attack.delta, defense.mask, pole_tol=pole_tol):
                attack_success[attack.label].add(defense.label)

    for defense in defenses:
        defense_success[defense.label] = set()
        for attack in attacks:
            if not is_destabilizing_static(A, B_w, C_r, attack.delta, defense.mask, pole_tol=pole_tol):
                defense_success[defense.label].add(attack.label)

    return StaticSuccessSets(
        attack_success=attack_success,
        defense_success=defense_success,
    )


def payoff_matrix_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    attacks: Sequence[mpm.AttackAction],
    defenses: Sequence[mpm.DefenseAction],
    pole_tol: float = 1e-8,
    ua: float = 1.0,
    ud: float = 0.0,
) -> np.ndarray:
    """Build the attacker payoff matrix for static destabilization checks."""
    U = np.zeros((len(attacks), len(defenses)), dtype=float)
    for i, attack in enumerate(attacks):
        for j, defense in enumerate(defenses):
            U[i, j] = (
                ua
                if is_destabilizing_static(A, B_w, C_r, attack.delta, defense.mask, pole_tol=pole_tol)
                else ud
            )
    return U


def eliminate_dominated_strategies_static(
    A: np.ndarray,
    B_w: np.ndarray,
    C_r: np.ndarray,
    attacks: Sequence[mpm.AttackAction],
    defenses: Sequence[mpm.DefenseAction],
    pole_tol: float = 1e-8,
) -> StaticReducedGame:
    """Iteratively remove strictly dominated actions."""
    cur_attacks = list(attacks)
    cur_defenses = list(defenses)

    changed = True
    while changed:
        changed = False
        success = compute_success_sets_static(A, B_w, C_r, cur_attacks, cur_defenses, pole_tol=pole_tol)
        dominated_attacks = mpm.dominated_attacks(success.attack_success)
        dominated_defenses = mpm.dominated_defenses(success.defense_success)

        if dominated_attacks:
            cur_attacks = [a for a in cur_attacks if a.label not in dominated_attacks]
            changed = True
        if dominated_defenses:
            cur_defenses = [d for d in cur_defenses if d.label not in dominated_defenses]
            changed = True

    payoff = payoff_matrix_static(A, B_w, C_r, cur_attacks, cur_defenses, pole_tol=pole_tol)
    return StaticReducedGame(payoff=payoff, attacks=cur_attacks, defenses=cur_defenses)


def defense_dataframe(
    defenses: Sequence[mpm.DefenseAction],
    c_w: Sequence[float],
    c_r: Sequence[float],
) -> pd.DataFrame:
    """Summarize rank-one defense masks."""
    rows: list[dict[str, float | int | str]] = []
    for defense in defenses:
        rows.append(
            {
                "defense": defense.label,
                "cost": mpm.defense_cost(mask=defense.mask, c_w=c_w, c_r=c_r),
                "num_defended_links": int((np.asarray(defense.mask) == 0).sum()),
                "num_uncovered_links": int((np.asarray(defense.mask) == 1).sum()),
            }
        )
    return pd.DataFrame(rows)


def payoff_dataframe(
    payoff: np.ndarray,
    attacks: Sequence[mpm.AttackAction],
    defenses: Sequence[mpm.DefenseAction],
) -> pd.DataFrame:
    """Return a labeled payoff matrix DataFrame."""
    return pd.DataFrame(
        np.asarray(payoff, dtype=float),
        index=[attack.label for attack in attacks],
        columns=[defense.label for defense in defenses],
    )


__all__ = [
    "HinfApproxResult",
    "StaticReducedGame",
    "StaticSuccessSets",
    "frequency_grid",
    "freqresp_matrix",
    "approximate_hinf_norm",
    "single_link_vulnerability_heatmap",
    "top_k_single_link_attacks",
    "rank1_static_delta_from_hinf",
    "single_link_static_delta",
    "masked_static_delta",
    "closed_loop_matrix_static",
    "closed_loop_eigs_static",
    "is_destabilizing_static",
    "compute_success_sets_static",
    "payoff_matrix_static",
    "eliminate_dominated_strategies_static",
    "defense_dataframe",
    "payoff_dataframe",
]
