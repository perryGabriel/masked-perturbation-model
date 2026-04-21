"""Chapter-2 feedback-vulnerability (FMP) utilities.

This module provides *numerical* evaluation helpers for small toy systems.
It does not provide symbolic/theorem proofs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal
import warnings

import numpy as np
import control

from .tf_tools import eye_tf, tf_matmul, tf_matrix

MatrixLike = np.ndarray
ThreatModel = Literal["full", "single_link"]
AccessModel = Literal["w0", "w1", "w2", "w3"]


@dataclass
class ContractSystem:
    """Small contract system for chapter-2 numerical exploration.

    Attributes
    ----------
    G:
        Contract matrix (typically static, small-dimensional for demos).
    alpha:
        Interconnection matrix with zero diagonal preferred.
    label:
        Optional identifier used in plots/reports.
    """

    G: MatrixLike
    alpha: MatrixLike
    label: str = "system"

    @property
    def n(self) -> int:
        return int(self.G.shape[0])


@dataclass
class VulnerabilityResult:
    value: float
    threat_model: ThreatModel
    access_model: str
    approximate: bool
    edge_warning: bool
    diagnostics: dict[str, float] | None = None


def build_contract_system(G_blocks: np.ndarray, alpha: np.ndarray, label: str = "system") -> ContractSystem:
    G = np.asarray(G_blocks, dtype=float)
    a = np.asarray(alpha, dtype=float)
    if G.ndim != 2 or G.shape[0] != G.shape[1]:
        raise ValueError("G must be square")
    if a.shape != G.shape:
        raise ValueError("alpha must have same shape as G")
    return ContractSystem(G=G, alpha=a, label=label)


def make_realization(G: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Compute P = (I - Q) G exactly for static matrices."""
    Gm = np.asarray(G, dtype=float)
    Qm = np.asarray(Q, dtype=float)
    I = np.eye(Gm.shape[0])
    return (I - Qm) @ Gm


def _is_dynamic_tf_matrix(x: object) -> bool:
    if isinstance(x, np.ndarray) and x.dtype == object:
        return any(isinstance(v, control.TransferFunction) for v in x.flat)
    return False


def make_realization_dynamic(G_tf: np.ndarray, Q_tf: np.ndarray) -> np.ndarray:
    """Compute ``P = (I - Q)G`` for transfer-matrix operands."""
    g_dyn = tf_matrix(G_tf)
    q_dyn = tf_matrix(Q_tf)
    if g_dyn.ndim != 2 or q_dyn.ndim != 2:
        raise TypeError("Dynamic realization expects 2D transfer-matrix operands for G_tf and Q_tf.")
    if g_dyn.shape[0] != g_dyn.shape[1]:
        raise ValueError(f"G_tf must be square; got shape {g_dyn.shape}.")
    if q_dyn.shape != g_dyn.shape:
        raise ValueError(f"Q_tf shape {q_dyn.shape} must match G_tf shape {g_dyn.shape}.")
    I = eye_tf(g_dyn.shape[0])
    i_minus_q = I.copy()
    for i in range(i_minus_q.shape[0]):
        for j in range(i_minus_q.shape[1]):
            i_minus_q[i, j] = i_minus_q[i, j] - q_dyn[i, j]
    return tf_matmul(i_minus_q, g_dyn)


def make_realization_checked(G: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Dispatch realization to static or dynamic path with explicit type checks."""
    g_is_dyn = _is_dynamic_tf_matrix(G)
    q_is_dyn = _is_dynamic_tf_matrix(Q)
    if g_is_dyn and q_is_dyn:
        return make_realization_dynamic(np.asarray(G, dtype=object), np.asarray(Q, dtype=object))
    if g_is_dyn != q_is_dyn:
        raise TypeError(
            "Cannot mix static and dynamic realization operands: "
            f"G is {'dynamic' if g_is_dyn else 'static'}, Q is {'dynamic' if q_is_dyn else 'static'}."
        )
    return make_realization(np.asarray(G, dtype=float), np.asarray(Q, dtype=float))


def compute_gamma(G: np.ndarray, alpha: np.ndarray, cond_max: float = 1e8) -> np.ndarray:
    """Compute Gamma = (I - G alpha)^(-1) for static matrices."""
    Gm = np.asarray(G, dtype=float)
    am = np.asarray(alpha, dtype=float)
    I = np.eye(Gm.shape[0])
    A = I - Gm @ am
    c = np.linalg.cond(A)
    if c > cond_max:
        raise ValueError(f"Ill-conditioned I - G alpha (cond={c:.2e})")
    return np.linalg.inv(A)


def compute_attack_map(Gamma: np.ndarray, Q: np.ndarray, Pbar: np.ndarray, Cbar: np.ndarray | None = None, cond_max: float = 1e8) -> np.ndarray:
    """Compute M_Q = Cbar * Gamma * (I - Q)^(-1) * Pbar.

    Exact for static matrix examples.
    """
    n = Gamma.shape[0]
    I = np.eye(n)
    A = I - np.asarray(Q, dtype=float)
    c = np.linalg.cond(A)
    if c > cond_max:
        raise ValueError(f"Ill-conditioned I - Q (cond={c:.2e})")
    inv_part = np.linalg.inv(A)
    C = np.eye(n) if Cbar is None else np.asarray(Cbar, dtype=float)
    return C @ np.asarray(Gamma, dtype=float) @ inv_part @ np.asarray(Pbar, dtype=float)


def _default_freq_grid() -> np.ndarray:
    # Broad default grid for continuous-time responses.
    return np.logspace(-3, 3, 400)


def _eval_tf_matrix_at_omega(M: np.ndarray, omega: float) -> np.ndarray:
    out = np.zeros(M.shape, dtype=complex)
    s = 1j * float(omega)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            out[i, j] = complex(control.evalfr(M[i, j], s))
    return out


def compute_attack_map_dynamic(
    Gamma: np.ndarray,
    Q_tf: np.ndarray,
    Pbar_tf: np.ndarray,
    Cbar: np.ndarray | None = None,
    freq_grid: np.ndarray | None = None,
    cond_max: float = 1e10,
) -> tuple[np.ndarray, dict[str, float]]:
    """Approximate dynamic attack-map samples on a frequency grid.

    Returns
    -------
    M_grid:
        Complex array with shape ``(n_out, n_in, n_freq)``.
    diagnostics:
        Includes grid coverage and worst conditioning of ``I - Q(jω)``.
    """
    q_dyn = tf_matrix(Q_tf)
    p_dyn = tf_matrix(Pbar_tf)
    if q_dyn.ndim != 2 or q_dyn.shape[0] != q_dyn.shape[1]:
        raise ValueError(f"Q_tf must be square 2D matrix; got shape {q_dyn.shape}.")
    if p_dyn.ndim != 2 or p_dyn.shape[0] != q_dyn.shape[0]:
        raise ValueError(
            f"Pbar_tf must have compatible shape (n, m) with n={q_dyn.shape[0]}; got {p_dyn.shape}."
        )
    omega = _default_freq_grid() if freq_grid is None else np.asarray(freq_grid, dtype=float).reshape(-1)
    if omega.size == 0:
        raise ValueError("freq_grid must contain at least one frequency sample.")
    C = np.eye(q_dyn.shape[0]) if Cbar is None else np.asarray(Cbar, dtype=float)
    Gm = np.asarray(Gamma, dtype=float)
    M_grid = np.zeros((C.shape[0], p_dyn.shape[1], omega.size), dtype=complex)
    cond_vals = np.zeros(omega.size, dtype=float)
    I = np.eye(q_dyn.shape[0])
    for k, w in enumerate(omega):
        q_w = _eval_tf_matrix_at_omega(q_dyn, w)
        p_w = _eval_tf_matrix_at_omega(p_dyn, w)
        A = I - q_w
        c = float(np.linalg.cond(A))
        cond_vals[k] = c
        if c > cond_max:
            raise ValueError(f"Ill-conditioned I - Q(jw) at w={w:.3e} (cond={c:.2e})")
        M_grid[:, :, k] = C @ Gm @ np.linalg.inv(A) @ p_w
    diagnostics = {
        "max_cond_i_minus_q": float(np.max(cond_vals)),
        "omega_min": float(np.min(omega)),
        "omega_max": float(np.max(omega)),
        "num_freq_samples": float(omega.size),
    }
    return M_grid, diagnostics


def access_matrix(system: ContractSystem, Q: np.ndarray, model: AccessModel) -> np.ndarray:
    """Construct Pbar for access models w0, w1, w2, w3."""
    n = system.n
    if model == "w0":
        return np.zeros((n, n))
    if model == "w1":
        return make_realization_checked(system.G, Q)
    if model == "w2":
        return np.eye(n)
    if model == "w3":
        return np.asarray(Q, dtype=float)
    raise ValueError(f"Unknown access model: {model}")


def well_posed(system: ContractSystem, Q: np.ndarray, cond_max: float = 1e8) -> bool:
    """Conservative static well-posedness/conditioning check."""
    try:
        _ = compute_gamma(system.G, system.alpha, cond_max=cond_max)
        _ = np.linalg.inv(np.eye(system.n) - np.asarray(Q, dtype=float))
    except Exception:
        return False
    return True


def vulnerability_full(system: ContractSystem, Q: np.ndarray, Pbar: np.ndarray, Cbar: np.ndarray | None = None, freq_grid: np.ndarray | None = None) -> VulnerabilityResult:
    """Measured vulnerability for full/unstructured attacker.

    Static matrices: exact spectral norm.
    Dynamic systems: would require grid approximation (not used in this static helper).
    """
    Gamma = compute_gamma(system.G, system.alpha)
    q_is_dyn = _is_dynamic_tf_matrix(Q)
    p_is_dyn = _is_dynamic_tf_matrix(Pbar)
    if q_is_dyn and p_is_dyn:
        M_grid, diagnostics = compute_attack_map_dynamic(Gamma, Q, Pbar, Cbar=Cbar, freq_grid=freq_grid)
        sigmas = np.array([np.linalg.svd(M_grid[:, :, k], compute_uv=False)[0] for k in range(M_grid.shape[2])], dtype=float)
        k_peak = int(np.argmax(sigmas))
        omega = _default_freq_grid() if freq_grid is None else np.asarray(freq_grid, dtype=float).reshape(-1)
        diagnostics["peak_frequency"] = float(omega[k_peak])
        diagnostics["grid_coverage_decades"] = float(np.log10(diagnostics["omega_max"] / diagnostics["omega_min"])) if diagnostics["omega_min"] > 0 else float("nan")
        return VulnerabilityResult(
            value=float(np.max(sigmas)),
            threat_model="full",
            access_model="direct",
            approximate=True,
            edge_warning=False,
            diagnostics=diagnostics,
        )
    if q_is_dyn != p_is_dyn:
        raise TypeError(
            "Cannot mix static and dynamic operands in vulnerability_full: "
            f"Q is {'dynamic' if q_is_dyn else 'static'}, Pbar is {'dynamic' if p_is_dyn else 'static'}."
        )
    M = compute_attack_map(Gamma, Q, Pbar, Cbar=Cbar)
    val = float(np.linalg.svd(M, compute_uv=False)[0])
    return VulnerabilityResult(value=val, threat_model="full", access_model="direct", approximate=False, edge_warning=False)


def vulnerability_single_link(system: ContractSystem, Q: np.ndarray, Pbar: np.ndarray, Cbar: np.ndarray | None = None, freq_grid: np.ndarray | None = None) -> VulnerabilityResult:
    """Measured vulnerability for single-link attacker.

    Static matrices: exact max-entry magnitude.
    """
    Gamma = compute_gamma(system.G, system.alpha)
    q_is_dyn = _is_dynamic_tf_matrix(Q)
    p_is_dyn = _is_dynamic_tf_matrix(Pbar)
    if q_is_dyn and p_is_dyn:
        M_grid, diagnostics = compute_attack_map_dynamic(Gamma, Q, Pbar, Cbar=Cbar, freq_grid=freq_grid)
        mags = np.array([np.max(np.abs(M_grid[:, :, k])) for k in range(M_grid.shape[2])], dtype=float)
        k_peak = int(np.argmax(mags))
        omega = _default_freq_grid() if freq_grid is None else np.asarray(freq_grid, dtype=float).reshape(-1)
        diagnostics["peak_frequency"] = float(omega[k_peak])
        diagnostics["grid_coverage_decades"] = float(np.log10(diagnostics["omega_max"] / diagnostics["omega_min"])) if diagnostics["omega_min"] > 0 else float("nan")
        return VulnerabilityResult(
            value=float(np.max(mags)),
            threat_model="single_link",
            access_model="direct",
            approximate=True,
            edge_warning=False,
            diagnostics=diagnostics,
        )
    if q_is_dyn != p_is_dyn:
        raise TypeError(
            "Cannot mix static and dynamic operands in vulnerability_single_link: "
            f"Q is {'dynamic' if q_is_dyn else 'static'}, Pbar is {'dynamic' if p_is_dyn else 'static'}."
        )
    M = compute_attack_map(Gamma, Q, Pbar, Cbar=Cbar)
    val = float(np.max(np.abs(M)))
    return VulnerabilityResult(value=val, threat_model="single_link", access_model="direct", approximate=False, edge_warning=False)


def q_strength_metric(Q: np.ndarray) -> float:
    """Proxy metric max_i sigma_max(Q_i), static simplification -> spectral norm."""
    return float(np.linalg.svd(np.asarray(Q, dtype=float), compute_uv=False)[0])
