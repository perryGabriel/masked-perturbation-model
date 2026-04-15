"""Projection-inspired design (finite-dimensional static surrogate).

This is a discretized/static surrogate inspired by chapter-2 projection ideas,
not an exact infinite-dimensional projection derivation.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .fmp import ContractSystem, access_matrix, compute_attack_map, compute_gamma, vulnerability_full, vulnerability_single_link, well_posed


@dataclass
class ProjectionIterationResult:
    q_vec: np.ndarray
    Q: np.ndarray
    surrogate_obj: float
    measured_vulnerability: float
    accepted: bool


@dataclass
class ProjectionDesignResult:
    Q_init: np.ndarray
    Q_final: np.ndarray
    iterations: list[ProjectionIterationResult]
    converged: bool
    rejected_steps: int


def _structure_indices(n: int, mask: np.ndarray | None = None) -> list[tuple[int, int]]:
    idx = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if mask is not None and not bool(mask[i, j]):
                continue
            idx.append((i, j))
    return idx


def _vec(M: np.ndarray) -> np.ndarray:
    return M.reshape(-1, order="F")


def projection_iteration(
    system: ContractSystem,
    Qk: np.ndarray,
    access_model: str = "w2",
    threat_model: str = "full",
    Pw: np.ndarray | None = None,
    structure_mask: np.ndarray | None = None,
    reg: float = 1e-4,
    damping: float = 0.5,
    cond_max: float = 1e6,
) -> ProjectionIterationResult:
    """One surrogate projection step based on regularized least squares."""
    n = system.n
    G = system.G
    Gamma = compute_gamma(G, system.alpha, cond_max=cond_max)
    I = np.eye(n)
    R = Gamma @ np.linalg.inv(I - Qk)
    Pbar = access_matrix(system, Qk, access_model)
    B = Pbar if Pw is None else (Pbar + Pw)

    idx = _structure_indices(n, structure_mask)
    A_cols = []
    for (i, j) in idx:
        E = np.zeros((n, n))
        E[i, j] = 1.0
        A_cols.append(_vec(R @ (-E @ G)))
    A = np.column_stack(A_cols) if A_cols else np.zeros((n * n, 0))
    b = _vec(R @ B)

    if A.shape[1] == 0:
        q_hat = np.array([])
    else:
        lhs = A.T @ A + reg * np.eye(A.shape[1])
        rhs = A.T @ b
        q_hat = np.linalg.solve(lhs, rhs)

    Q_hat = np.zeros_like(Qk)
    for k, (i, j) in enumerate(idx):
        Q_hat[i, j] = q_hat[k]

    Q_new = (1 - damping) * Qk + damping * Q_hat
    accepted = well_posed(system, Q_new, cond_max=cond_max)
    if not accepted:
        Q_new = Qk.copy()

    surrogate_obj = float(np.linalg.norm(R @ (-Q_new @ G + B), ord="fro"))

    if threat_model == "full":
        v = vulnerability_full(system, Q_new, access_matrix(system, Q_new, access_model)).value
    elif threat_model == "single_link":
        v = vulnerability_single_link(system, Q_new, access_matrix(system, Q_new, access_model)).value
    else:
        raise ValueError("threat_model must be 'full' or 'single_link'")

    return ProjectionIterationResult(
        q_vec=q_hat,
        Q=Q_new,
        surrogate_obj=surrogate_obj,
        measured_vulnerability=float(v),
        accepted=accepted,
    )


def projection_design(
    system: ContractSystem,
    Q0: np.ndarray | None = None,
    max_iter: int = 20,
    tol: float = 1e-6,
    access_model: str = "w2",
    threat_model: str = "full",
    structure_mask: np.ndarray | None = None,
    reg: float = 1e-4,
    damping: float = 0.5,
    cond_max: float = 1e6,
) -> ProjectionDesignResult:
    n = system.n
    Qk = np.zeros((n, n)) if Q0 is None else np.asarray(Q0, dtype=float).copy()
    history: list[ProjectionIterationResult] = []
    rejected = 0

    for _ in range(max_iter):
        step = projection_iteration(
            system=system,
            Qk=Qk,
            access_model=access_model,
            threat_model=threat_model,
            structure_mask=structure_mask,
            reg=reg,
            damping=damping,
            cond_max=cond_max,
        )
        history.append(step)
        if not step.accepted:
            rejected += 1
        if np.linalg.norm(step.Q - Qk, ord="fro") < tol:
            return ProjectionDesignResult(Q_init=np.zeros((n, n)) if Q0 is None else np.asarray(Q0, dtype=float), Q_final=step.Q, iterations=history, converged=True, rejected_steps=rejected)
        Qk = step.Q

    return ProjectionDesignResult(Q_init=np.zeros((n, n)) if Q0 is None else np.asarray(Q0, dtype=float), Q_final=Qk, iterations=history, converged=False, rejected_steps=rejected)
