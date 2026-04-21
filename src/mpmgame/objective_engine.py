from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from .fmp import ContractSystem, access_matrix, vulnerability_full, vulnerability_single_link
from .dsf_param import DynamicQ

ThreatModel = Literal["full", "single_link"]
ObjectiveType = Literal["true", "frobenius_proxy", "lp_linearized_proxy"]


@dataclass
class ParameterizationSpec:
    """Finite-dimensional map between theta and Q for benchmark runs."""

    kind: str = "static_hollow"
    mask: np.ndarray | None = None
    bounds: tuple[float, float] | None = None
    g_hat: float | None = None
    q_num_degree: int = 0
    q_den_degree: int = 0
    shared_denominator: bool = True
    stability_parameterization: str = "fixed_den"
    zero_diagonal: bool = True
    freq_grid: np.ndarray | None = None
    freq_weighting: np.ndarray | None = None


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
    spec = problem.parameterization
    n_free = len(free_entries(problem.system.n, spec.mask))
    if spec.kind == "static_hollow":
        return n_free
    if spec.kind == "dsf_poly":
        num_params = n_free * (int(spec.q_num_degree) + 1)
        den_degree = int(spec.q_den_degree)
        if den_degree > 0 and spec.stability_parameterization != "fixed_den":
            if spec.shared_denominator:
                den_params = den_degree
            else:
                den_params = n_free * den_degree
        else:
            den_params = 0
        return int(num_params + den_params)
    raise ValueError(f"Unsupported parameterization.kind={spec.kind}")


def _dsf_poly_layout(problem: ProblemSpec) -> tuple[list[tuple[int, int]], int, int]:
    spec = problem.parameterization
    idx = free_entries(problem.system.n, spec.mask)
    return idx, int(spec.q_num_degree) + 1, int(spec.q_den_degree)


def theta_to_q(theta: np.ndarray, problem: ProblemSpec) -> tuple[np.ndarray, dict[str, Any]]:
    n = problem.system.n
    spec = problem.parameterization
    idx = free_entries(n, spec.mask)
    t = np.asarray(theta, dtype=float).reshape(-1)
    if spec.kind == "static_hollow":
        if len(t) != len(idx):
            raise ValueError(f"theta length {len(t)} != required {len(idx)}")
        Q = np.zeros((n, n), dtype=float)
        for k, (i, j) in enumerate(idx):
            Q[i, j] = t[k]
        return Q, {"n_free": len(idx), "indices": idx}

    if spec.kind == "dsf_poly":
        idx, num_terms, den_degree = _dsf_poly_layout(problem)
        n_free = len(idx)
        num_count = n_free * num_terms
        den_count = 0
        if den_degree > 0 and spec.stability_parameterization != "fixed_den":
            den_count = den_degree if spec.shared_denominator else n_free * den_degree
        required = num_count + den_count
        if len(t) != required:
            raise ValueError(f"theta length {len(t)} != required {required}")

        Q = np.zeros((n, n), dtype=float)
        num_coeffs = np.zeros((n, n, num_terms), dtype=float)
        cursor = 0
        for i, j in idx:
            coeffs = t[cursor : cursor + num_terms]
            num_coeffs[i, j, :] = coeffs
            Q[i, j] = coeffs[0]
            cursor += num_terms

        den_coeffs = None
        if den_count > 0:
            if spec.shared_denominator:
                den_coeffs = t[cursor : cursor + den_degree].copy()
            else:
                den_coeffs = np.zeros((n, n, den_degree), dtype=float)
                for i, j in idx:
                    den_coeffs[i, j, :] = t[cursor : cursor + den_degree]
                    cursor += den_degree

        if spec.zero_diagonal:
            np.fill_diagonal(Q, 0.0)
            for d in range(num_terms):
                np.fill_diagonal(num_coeffs[:, :, d], 0.0)

        return Q, {
            "n_free": n_free,
            "indices": idx,
            "kind": "dsf_poly",
            "q_num_degree": spec.q_num_degree,
            "q_den_degree": spec.q_den_degree,
            "num_coeffs": num_coeffs,
            "den_coeffs": den_coeffs,
        }

    raise ValueError(f"Unsupported parameterization.kind={spec.kind}")


def q_to_theta(Q: np.ndarray, problem: ProblemSpec) -> np.ndarray:
    q = np.asarray(Q, dtype=float)
    spec = problem.parameterization
    idx = free_entries(problem.system.n, spec.mask)
    if spec.kind == "static_hollow":
        return np.array([q[i, j] for i, j in idx], dtype=float)
    if spec.kind == "dsf_poly":
        _, num_terms, den_degree = _dsf_poly_layout(problem)
        vals: list[float] = []
        for i, j in idx:
            vals.append(float(q[i, j]))
            if num_terms > 1:
                vals.extend([0.0] * (num_terms - 1))
        if den_degree > 0 and spec.stability_parameterization != "fixed_den":
            den_count = den_degree if spec.shared_denominator else len(idx) * den_degree
            vals.extend([0.0] * den_count)
        return np.array(vals, dtype=float)
    raise ValueError(f"Unsupported parameterization.kind={spec.kind}")


def _effective_degree(coeffs: np.ndarray, tol: float = 1e-12) -> int:
    c = np.asarray(coeffs, dtype=float).reshape(-1)
    nz = np.flatnonzero(np.abs(c) > tol)
    return int(nz[-1]) if nz.size else 0


def _eval_dsf_ratio(num_asc: np.ndarray, den_tail_asc: np.ndarray | None, omega: float) -> complex:
    z_inv = np.exp(-1j * omega)
    num = np.sum(np.asarray(num_asc, dtype=float) * (z_inv ** np.arange(len(num_asc))))
    den = 1.0
    if den_tail_asc is not None and len(den_tail_asc) > 0:
        den = den + np.sum(np.asarray(den_tail_asc, dtype=float) * (z_inv ** np.arange(1, len(den_tail_asc) + 1)))
    return complex(num / den)


def _check_feasibility(
    Q: np.ndarray, problem: ProblemSpec, meta: dict[str, Any] | None = None
) -> tuple[bool, list[str], dict[str, Any]]:
    reasons: list[str] = []
    n = problem.system.n
    cond_val: float | None = None
    cond_vals: list[float] = []

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

    if problem.parameterization.kind == "dsf_poly":
        if meta is None or "num_coeffs" not in meta:
            reasons.append("missing_dsf_metadata")
        else:
            spec = problem.parameterization
            idx, _, _ = _dsf_poly_layout(problem)
            num_coeffs = np.asarray(meta["num_coeffs"], dtype=float)
            den_coeffs = None if meta.get("den_coeffs", None) is None else np.asarray(meta["den_coeffs"], dtype=float)
            dyn_q = DynamicQ.from_tensors(
                num_coeffs=num_coeffs,
                den_coeffs=den_coeffs,
                mask=spec.mask,
                zero_diagonal=spec.zero_diagonal,
                metadata={"problem_id": problem.problem_id},
                enforce=True,
            )
            reasons.extend(dyn_q.verify_constraints())

            den_is_optimized = int(spec.q_den_degree) > 0 and spec.stability_parameterization != "fixed_den"
            if int(spec.q_den_degree) >= 0 and int(spec.q_num_degree) > int(spec.q_den_degree):
                reasons.append("improper_entry")
            else:
                for i, j in idx:
                    num_deg = _effective_degree(num_coeffs[i, j, :])
                    if den_coeffs is None:
                        den_deg = 0
                    elif den_coeffs.ndim == 1:
                        den_deg = _effective_degree(np.concatenate([[1.0], den_coeffs]))
                    else:
                        den_deg = _effective_degree(np.concatenate([[1.0], den_coeffs[i, j, :]]))
                    if num_deg > den_deg:
                        reasons.append("improper_entry")
                        break

            if den_is_optimized and den_coeffs is not None:
                if den_coeffs.ndim == 1:
                    roots = np.roots(np.concatenate([[1.0], den_coeffs]))
                    if np.any(np.abs(roots) >= 1.0 - 1e-10):
                        reasons.append("unstable_denominator")
                else:
                    unstable = False
                    for i, j in idx:
                        roots = np.roots(np.concatenate([[1.0], den_coeffs[i, j, :]]))
                        if np.any(np.abs(roots) >= 1.0 - 1e-10):
                            unstable = True
                            break
                    if unstable:
                        reasons.append("unstable_denominator")

            freq_grid = problem.freq_grid
            if freq_grid is None:
                freq_grid = problem.parameterization.freq_grid

            if freq_grid is not None and len(freq_grid) > 0:
                for w in np.asarray(freq_grid, dtype=float).reshape(-1):
                    try:
                        Qw = np.zeros((n, n), dtype=complex)
                        for i, j in idx:
                            den_tail = None if den_coeffs is None else (den_coeffs if den_coeffs.ndim == 1 else den_coeffs[i, j, :])
                            Qw[i, j] = _eval_dsf_ratio(num_coeffs[i, j, :], den_tail, float(w))
                        A = np.eye(n, dtype=complex) - Qw
                        c = float(np.linalg.cond(A))
                        cond_vals.append(c)
                    except Exception:
                        reasons.append("freq_inversion_failure")
                        break
                if cond_vals:
                    cond_val = float(np.max(cond_vals))
                    if cond_val > problem.cond_max:
                        reasons.append("freq_ill_conditioned")
            else:
                try:
                    A = np.eye(n) - Q
                    cond_val = float(np.linalg.cond(A))
                    cond_vals.append(cond_val)
                    if cond_val > problem.cond_max:
                        reasons.append("ill_conditioned")
                    _ = np.linalg.inv(A)
                except Exception:
                    reasons.append("inversion_failure")
    else:
        try:
            A = np.eye(n) - Q
            cond_val = float(np.linalg.cond(A))
            cond_vals.append(cond_val)
            if cond_val > problem.cond_max:
                reasons.append("ill_conditioned")
            _ = np.linalg.inv(A)
        except Exception:
            reasons.append("inversion_failure")

    dedup_reasons = list(dict.fromkeys(reasons))
    return len(dedup_reasons) == 0, dedup_reasons, {"cond_i_minus_q": cond_val, "cond_i_minus_q_grid": cond_vals}


def evaluate_theta(theta: np.ndarray, problem: ProblemSpec) -> EvalResult:
    Q, meta = theta_to_q(theta, problem)
    feasible, reasons, diag = _check_feasibility(Q, problem, meta=meta)

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
        q_for_access: np.ndarray = Q
        if (
            problem.parameterization.kind == "dsf_poly"
            and (int(problem.parameterization.q_num_degree) > 0 or int(problem.parameterization.q_den_degree) > 0)
            and "num_coeffs" in meta
        ):
            q_for_access = DynamicQ.from_tensors(
                num_coeffs=np.asarray(meta["num_coeffs"], dtype=float),
                den_coeffs=None if meta.get("den_coeffs", None) is None else np.asarray(meta["den_coeffs"], dtype=float),
                mask=problem.parameterization.mask,
                zero_diagonal=problem.parameterization.zero_diagonal,
                metadata={"problem_id": problem.problem_id},
            ).to_tf_matrix()

        pbar = access_matrix(problem.system, q_for_access, problem.access_model)
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
