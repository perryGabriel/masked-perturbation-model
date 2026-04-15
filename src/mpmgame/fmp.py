"""Chapter-2 feedback-vulnerability (FMP) utilities.

This module provides *numerical* evaluation helpers for small toy systems.
It does not provide symbolic/theorem proofs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal
import warnings

import numpy as np

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


def access_matrix(system: ContractSystem, Q: np.ndarray, model: AccessModel) -> np.ndarray:
    """Construct Pbar for access models w0, w1, w2, w3."""
    n = system.n
    if model == "w0":
        return np.zeros((n, n))
    if model == "w1":
        return make_realization(system.G, Q)
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
    M = compute_attack_map(Gamma, Q, Pbar, Cbar=Cbar)
    val = float(np.linalg.svd(M, compute_uv=False)[0])
    return VulnerabilityResult(value=val, threat_model="full", access_model="direct", approximate=False, edge_warning=False)


def vulnerability_single_link(system: ContractSystem, Q: np.ndarray, Pbar: np.ndarray, Cbar: np.ndarray | None = None, freq_grid: np.ndarray | None = None) -> VulnerabilityResult:
    """Measured vulnerability for single-link attacker.

    Static matrices: exact max-entry magnitude.
    """
    Gamma = compute_gamma(system.G, system.alpha)
    M = compute_attack_map(Gamma, Q, Pbar, Cbar=Cbar)
    val = float(np.max(np.abs(M)))
    return VulnerabilityResult(value=val, threat_model="single_link", access_model="direct", approximate=False, edge_warning=False)


def q_strength_metric(Q: np.ndarray) -> float:
    """Proxy metric max_i sigma_max(Q_i), static simplification -> spectral norm."""
    return float(np.linalg.svd(np.asarray(Q, dtype=float), compute_uv=False)[0])
