from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from .fmp import ContractSystem, access_matrix, vulnerability_full, vulnerability_single_link

ThreatModel = Literal["full", "single_link"]
ObjectiveType = Literal["true", "frobenius_proxy", "lp_linearized_proxy"]


@dataclass
class ParameterizationSpec:
    """Finite-dimensional map between theta and Q for benchmark runs."""

    kind: str = "static_hollow"
    mask: np.ndarray | None = None
    bounds: tuple[float, float] | None = None
    g_hat: float | None = None


@dataclass
class ProblemSpec:
    problem_id: str
    system: ContractSystem
    access_model: str = "w2"
    threat_model: ThreatModel = "full"
    parameterization: ParameterizationSpec = field(default_factory=ParameterizationSpec)
    objective_type: ObjectiveType = "true"
    freq_grid: np.ndarray | None = None
    cond_max: float = 1e8
    infeasible_penalty: float = 1e6
    hard_reject: bool = False


@dataclass
class EvalResult:
    objective: float
    true_objective: float
    surrogate_objective: float | None
    feasible: bool
    penalty: float
    Q: np.ndarray
    diagnostics: dict[str, Any]


@dataclass
class OptimizationResult:
    algorithm: str
    problem_id: str
    restart_id: int
    seed: int
    initialization: str
    success: bool
    feasible: bool
    timeout: bool
    objective: float
    true_objective: float
    surrogate_objective: float | None
    runtime_sec: float
    nfev: int | None
    nit: int | None
    param_norm: float
    cond_i_minus_q: float | None
    message: str
    theta: np.ndarray
    Q: np.ndarray
    diagnostics: dict[str, Any] = field(default_factory=dict)
    trace: list[dict[str, float]] = field(default_factory=list)


def free_entries(n: int, mask: np.ndarray | None = None) -> list[tuple[int, int]]:
    idx: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if mask is not None and not bool(mask[i, j]):
                continue
            idx.append((i, j))
    return idx


def theta_dim(problem: ProblemSpec) -> int:
    return len(free_entries(problem.system.n, problem.parameterization.mask))


def theta_to_q(theta: np.ndarray, problem: ProblemSpec) -> tuple[np.ndarray, dict[str, Any]]:
    n = problem.system.n
    idx = free_entries(n, problem.parameterization.mask)
    t = np.asarray(theta, dtype=float).reshape(-1)
    if len(t) != len(idx):
        raise ValueError(f"theta length {len(t)} != required {len(idx)}")
    Q = np.zeros((n, n), dtype=float)
    for k, (i, j) in enumerate(idx):
        Q[i, j] = t[k]
    return Q, {"n_free": len(idx), "indices": idx}


def q_to_theta(Q: np.ndarray, problem: ProblemSpec) -> np.ndarray:
    q = np.asarray(Q, dtype=float)
    idx = free_entries(problem.system.n, problem.parameterization.mask)
    return np.array([q[i, j] for i, j in idx], dtype=float)


def _check_feasibility(Q: np.ndarray, problem: ProblemSpec) -> tuple[bool, list[str], dict[str, Any]]:
    reasons: list[str] = []
    n = problem.system.n
    I = np.eye(n)
    cond_val: float | None = None

    if not np.allclose(np.diag(Q), 0.0):
        reasons.append("non_hollow")

    if problem.parameterization.mask is not None:
        forbidden = (~problem.parameterization.mask.astype(bool)) & (~np.eye(n, dtype=bool))
        if np.any(np.abs(Q[forbidden]) > 1e-10):
            reasons.append("mask_violation")

    if problem.parameterization.g_hat is not None:
        qnorm = float(np.linalg.norm(Q, ord=2))
        if qnorm > float(problem.parameterization.g_hat) + 1e-12:
            reasons.append("gain_bound_violation")

    bounds = problem.parameterization.bounds
    if bounds is not None:
        lo, hi = bounds
        offdiag = Q[~np.eye(n, dtype=bool)]
        if np.any(offdiag < lo - 1e-12) or np.any(offdiag > hi + 1e-12):
            reasons.append("box_bound_violation")

    try:
        A = I - Q
        cond_val = float(np.linalg.cond(A))
        if cond_val > problem.cond_max:
            reasons.append("ill_conditioned")
        _ = np.linalg.inv(A)
    except Exception:
        reasons.append("inversion_failure")

    return len(reasons) == 0, reasons, {"cond_i_minus_q": cond_val}


def evaluate_theta(theta: np.ndarray, problem: ProblemSpec) -> EvalResult:
    Q, meta = theta_to_q(theta, problem)
    feasible, reasons, diag = _check_feasibility(Q, problem)

    penalty = 0.0
    if not feasible:
        penalty = problem.infeasible_penalty * (1.0 + 0.1 * len(reasons))
        if problem.hard_reject:
            return EvalResult(
                objective=penalty,
                true_objective=np.inf,
                surrogate_objective=None,
                feasible=False,
                penalty=penalty,
                Q=Q,
                diagnostics={**meta, **diag, "reasons": reasons},
            )

    true_obj = np.inf
    surrogate_obj: float | None = None
    try:
        pbar = access_matrix(problem.system, Q, problem.access_model)
        if problem.threat_model == "full":
            true_obj = vulnerability_full(problem.system, Q, pbar, freq_grid=problem.freq_grid).value
        else:
            true_obj = vulnerability_single_link(problem.system, Q, pbar, freq_grid=problem.freq_grid).value

        if problem.objective_type == "true":
            obj = true_obj
        elif problem.objective_type == "frobenius_proxy":
            M = np.linalg.inv(np.eye(problem.system.n) - Q) @ pbar
            surrogate_obj = float(np.linalg.norm(M, ord="fro"))
            obj = surrogate_obj
        elif problem.objective_type == "lp_linearized_proxy":
            M_lin = pbar + Q @ pbar
            surrogate_obj = float(np.max(np.abs(M_lin)))
            obj = surrogate_obj
        else:
            raise ValueError(f"Unsupported objective_type={problem.objective_type}")
    except Exception as exc:
        feasible = False
        reasons.append(f"objective_failure:{type(exc).__name__}")
        obj = problem.infeasible_penalty
        true_obj = np.inf

    if not np.isfinite(obj):
        feasible = False
        reasons.append("nan_inf_objective")
        obj = problem.infeasible_penalty

    if not feasible:
        obj = float(obj + penalty)

    return EvalResult(
        objective=float(obj),
        true_objective=float(true_obj),
        surrogate_objective=surrogate_obj,
        feasible=feasible,
        penalty=float(penalty),
        Q=Q,
        diagnostics={**meta, **diag, "reasons": reasons},
    )
