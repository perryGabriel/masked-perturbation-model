from __future__ import annotations

import dataclasses
import json
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
try:
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else []

from .design_lp import lp_relaxation_design
from .design_projection import projection_design
from .fmp import build_contract_system
from .initialization import generate_initializations
from .objective_engine import OptimizationResult, ParameterizationSpec, ProblemSpec, evaluate_theta, q_to_theta, theta_dim
from .optimize_nonlinear import OPTIONAL_ALGORITHMS, REQUIRED_ALGORITHMS, available_optional_algorithms, optimize_with_algorithm
from .reporting_benchmarks import (
    aggregate_results,
    generate_plots,
    write_incremental_result_markdown,
    write_report_markdown,
)
from .timeouts import run_with_timeout


def _toy_problem(problem_id: str, n: int, access_model: str, threat_model: str) -> ProblemSpec:
    G = np.diag(np.linspace(0.4, 0.7, n))
    alpha = np.zeros((n, n))
    for i in range(n - 1):
        alpha[i, i + 1] = 0.1
        alpha[i + 1, i] = 0.05
    system = build_contract_system(G, alpha, label=problem_id)
    param = ParameterizationSpec(kind="static_hollow", bounds=(-0.25, 0.25), g_hat=0.6)
    return ProblemSpec(problem_id=problem_id, system=system, access_model=access_model, threat_model=threat_model, parameterization=param)


def _toy_problem_dsf(
    problem_id: str,
    n: int,
    access_model: str,
    threat_model: str,
    *,
    q_num_degree: int = 1,
    q_den_degree: int = 1,
    shared_denominator: bool = True,
    stability_parameterization: str = "schur_reflect",
) -> ProblemSpec:
    G = np.diag(np.linspace(0.4, 0.7, n))
    alpha = np.zeros((n, n))
    for i in range(n - 1):
        alpha[i, i + 1] = 0.1
        alpha[i + 1, i] = 0.05
    system = build_contract_system(G, alpha, label=problem_id)
    param = ParameterizationSpec(
        kind="dsf_poly",
        bounds=(-0.25, 0.25),
        g_hat=0.6,
        q_num_degree=q_num_degree,
        q_den_degree=q_den_degree,
        shared_denominator=shared_denominator,
        stability_parameterization=stability_parameterization,
    )
    return ProblemSpec(problem_id=problem_id, system=system, access_model=access_model, threat_model=threat_model, parameterization=param)


def benchmark_problem_registry(preset: str = "quick") -> list[ProblemSpec]:
    quick = [
        _toy_problem("toy2_w1_full", 2, "w1", "full"),
        _toy_problem("toy3_w2_single", 3, "w2", "single_link"),
        _toy_problem("toy3_w3_full", 3, "w3", "full"),
    ]
    if preset == "quick":
        return quick
    return quick + [
        _toy_problem("toy2_w2_single", 2, "w2", "single_link"),
        _toy_problem("toy3_w1_full", 3, "w1", "full"),
    ]


def benchmark_problem_registry_dsf(preset: str = "quick") -> list[ProblemSpec]:
    quick = [
        _toy_problem_dsf("toy2_w1_full_dsf", 2, "w1", "full"),
        _toy_problem_dsf("toy3_w2_single_dsf", 3, "w2", "single_link"),
        _toy_problem_dsf("toy3_w3_full_dsf", 3, "w3", "full"),
    ]
    if preset == "quick":
        return quick
    return quick + [
        _toy_problem_dsf("toy2_w2_single_dsf", 2, "w2", "single_link"),
        _toy_problem_dsf("toy3_w1_full_dsf", 3, "w1", "full"),
    ]


def _to_row(r: OptimizationResult) -> dict:
    d = dataclasses.asdict(r)
    diagnostics = dict(r.diagnostics or {})
    param_kind = str(diagnostics.get("kind", "static_hollow"))
    if param_kind == "":
        param_kind = "static_hollow"
    d["theta"] = np.asarray(r.theta).tolist()
    d["Q"] = np.asarray(r.Q).tolist()
    d["param_kind"] = param_kind
    d["Q_coeffs"] = None
    d["Q_meta"] = None
    d["P_coeffs"] = None

    if param_kind == "dsf_poly":
        num_coeffs = diagnostics.get("num_coeffs", None)
        den_coeffs = diagnostics.get("den_coeffs", None)
        q_payload: dict[str, object] = {"kind": "dsf_poly", "num_coeffs": None, "den_coeffs": None}
        if num_coeffs is not None:
            q_payload["num_coeffs"] = np.asarray(num_coeffs, dtype=float).tolist()
        if den_coeffs is not None:
            q_payload["den_coeffs"] = np.asarray(den_coeffs, dtype=float).tolist()
        d["Q_coeffs"] = q_payload

        q_meta = {
            "q_num_degree": diagnostics.get("q_num_degree"),
            "q_den_degree": diagnostics.get("q_den_degree"),
            "n_free": diagnostics.get("n_free"),
            "indices": diagnostics.get("indices"),
            "cond_i_minus_q_grid": diagnostics.get("cond_i_minus_q_grid"),
            "reasons": diagnostics.get("reasons"),
        }
        d["Q_meta"] = {k: v for k, v in q_meta.items() if v is not None}

        p_coeffs = diagnostics.get("P_coeffs", None)
        if p_coeffs is not None:
            d["P_coeffs"] = p_coeffs
    return d


def _baseline_results(problem: ProblemSpec) -> list[OptimizationResult]:
    rows: list[OptimizationResult] = []
    zero = np.zeros(theta_dim(problem))
    z = evaluate_theta(zero, problem)
    rows.append(
        OptimizationResult(
            algorithm="baseline_zero", problem_id=problem.problem_id, restart_id=0, seed=0, initialization="zero",
            success=True, feasible=z.feasible, timeout=False, objective=z.objective, true_objective=z.true_objective,
            surrogate_objective=z.surrogate_objective, runtime_sec=0.0, nfev=1, nit=0, param_norm=float(np.linalg.norm(zero)),
            cond_i_minus_q=z.diagnostics.get("cond_i_minus_q"), message="fixed baseline", theta=zero, Q=z.Q, diagnostics=z.diagnostics
        )
    )

    t0 = time.perf_counter()
    try:
        pres = projection_design(problem.system, threat_model=problem.threat_model, access_model=problem.access_model, structure_mask=problem.parameterization.mask)
        ptheta = q_to_theta(pres.Q_final, problem)
        pev = evaluate_theta(ptheta, problem)
        msg = f"converged={pres.converged}"
        success = True
    except Exception as exc:
        ptheta = zero.copy(); pev = evaluate_theta(ptheta, problem); msg = f"error:{exc}"; success = False
    rows.append(OptimizationResult("baseline_projection", problem.problem_id, 0, 0, "projection", success, pev.feasible, False, pev.objective, pev.true_objective, pev.surrogate_objective, time.perf_counter()-t0, None, None, float(np.linalg.norm(ptheta)), pev.diagnostics.get("cond_i_minus_q"), msg, ptheta, pev.Q, diagnostics=pev.diagnostics))

    t1 = time.perf_counter()
    try:
        lres = lp_relaxation_design(problem.system, access_model=problem.access_model, threat_model=problem.threat_model, g_hat=problem.parameterization.g_hat or 0.2, structure_mask=problem.parameterization.mask)
        ltheta = q_to_theta(lres.Q, problem)
        lev = evaluate_theta(ltheta, problem)
        msg = lres.message
        success = lres.success
    except Exception as exc:
        ltheta = zero.copy(); lev = evaluate_theta(ltheta, problem); msg = f"error:{exc}"; success = False
    rows.append(OptimizationResult("baseline_lp", problem.problem_id, 0, 0, "lp", success, lev.feasible, False, lev.objective, lev.true_objective, lev.surrogate_objective, time.perf_counter()-t1, None, None, float(np.linalg.norm(ltheta)), lev.diagnostics.get("cond_i_minus_q"), msg, ltheta, lev.Q, diagnostics=lev.diagnostics))
    return rows


def _execute_run(problem: ProblemSpec, algorithm: str, init_theta: np.ndarray, seed: int, restart_id: int, init_name: str, maxiter: int, maxfun: int):
    return optimize_with_algorithm(problem, algorithm, init_theta, seed=seed, restart_id=restart_id, initialization=init_name, maxiter=maxiter, maxfun=maxfun)


def run_benchmark_suite(
    problem_specs: Iterable[ProblemSpec],
    algorithms: list[str] | None = None,
    n_restarts: int = 5,
    timeout_per_run_sec: float = 45.0,
    algorithm_timeout_sec: float | None = None,
    output_dir: str | Path = "results/optimizer_benchmarks",
    report_dir: str | Path = "reports/optimizer_benchmarks",
    seed: int = 7,
    show_progress: bool = True,
    maxiter: int = 200,
    maxfun: int = 2000,
    incremental_notes_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    algos = algorithms or ["Nelder-Mead", "Powell", "SLSQP", "L-BFGS-B", "differential_evolution"]
    optional_avail = available_optional_algorithms()
    rows: list[OptimizationResult] = []

    probs = list(problem_specs)
    pbar_prob = tqdm(probs, disable=not show_progress, desc="Problems")
    for prob in pbar_prob:
        rows.extend(_baseline_results(prob))

        pbar_algo = tqdm(algos, disable=not show_progress, leave=False, desc=f"Algorithms:{prob.problem_id}")
        for algo in pbar_algo:
            algo_start = time.perf_counter()
            if algo in OPTIONAL_ALGORITHMS and not optional_avail.get(algo, False):
                rows.append(OptimizationResult(algo, prob.problem_id, -1, seed, "n/a", False, False, False, np.inf, np.inf, None, 0.0, 0, 0, 0.0, None, "optional dependency unavailable", np.zeros(theta_dim(prob)), np.zeros((prob.system.n, prob.system.n))))
                continue

            inits = generate_initializations(prob, n_random=max(1, n_restarts // 2), seed=seed)
            restart_points = inits[:n_restarts]
            feasible_count = 0
            timeout_count = 0
            current_best = np.inf
            pbar_restart = tqdm(enumerate(restart_points), total=len(restart_points), disable=not show_progress, leave=False, desc=f"Restarts:{algo}")
            for restart_id, init in pbar_restart:
                run_seed = seed + restart_id * 13
                out = run_with_timeout(_execute_run, timeout_per_run_sec, prob, algo, init.theta, run_seed, restart_id, init.name, maxiter, maxfun)
                if out.timed_out:
                    timeout_count += 1
                    rows.append(OptimizationResult(algo, prob.problem_id, restart_id, run_seed, init.name, False, False, True, np.inf, np.inf, None, timeout_per_run_sec, None, None, float(np.linalg.norm(init.theta)), None, "timeout", init.theta, np.zeros((prob.system.n, prob.system.n))))
                elif out.error is not None:
                    rows.append(OptimizationResult(algo, prob.problem_id, restart_id, run_seed, init.name, False, False, False, np.inf, np.inf, None, 0.0, None, None, float(np.linalg.norm(init.theta)), None, out.error, init.theta, np.zeros((prob.system.n, prob.system.n))))
                else:
                    res: OptimizationResult = out.value
                    rows.append(res)
                    if res.feasible:
                        feasible_count += 1
                        current_best = min(current_best, res.true_objective)
                if incremental_notes_dir is not None:
                    partial_df = pd.DataFrame([_to_row(r) for r in rows])
                    note_path = Path(incremental_notes_dir) / prob.problem_id / f"{algo}.md"
                    write_incremental_result_markdown(note_path, partial_df, prob.problem_id, algo)
                if hasattr(pbar_restart, "set_postfix"):
                    pbar_restart.set_postfix(best=f"{current_best:.3g}" if np.isfinite(current_best) else "inf", feasible=feasible_count, timeout=timeout_count)

                if algorithm_timeout_sec is not None and (time.perf_counter() - algo_start) > algorithm_timeout_sec:
                    break

    raw_df = pd.DataFrame([_to_row(r) for r in rows])
    summary_df = aggregate_results(raw_df)

    out_path = Path(output_dir); out_path.mkdir(parents=True, exist_ok=True)
    rep_path = Path(report_dir); rep_path.mkdir(parents=True, exist_ok=True)
    raw_csv = out_path / "benchmark_raw_results.csv"
    summary_csv = out_path / "benchmark_summary.csv"
    raw_df.to_csv(raw_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)

    generate_plots(raw_df, summary_df, rep_path)
    write_report_markdown(rep_path / "benchmark_report.md", raw_df, summary_df, setup_text=f"algorithms={algos}, n_restarts={n_restarts}, timeout_per_run_sec={timeout_per_run_sec}")
    (out_path / "run_config.json").write_text(json.dumps({"algorithms": algos, "n_restarts": n_restarts, "timeout_per_run_sec": timeout_per_run_sec, "seed": seed}, indent=2))
    return raw_df, summary_df
