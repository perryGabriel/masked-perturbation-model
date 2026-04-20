from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.optimize import basinhopping, differential_evolution, dual_annealing, minimize, shgo

from .initialization import clip_theta_to_bounds
from .objective_engine import OptimizationResult, ProblemSpec, evaluate_theta

REQUIRED_ALGORITHMS = [
    "Nelder-Mead",
    "Powell",
    "COBYLA",
    "SLSQP",
    "L-BFGS-B",
    "trust-constr",
    "differential_evolution",
    "dual_annealing",
    "basinhopping",
    "shgo",
]

OPTIONAL_ALGORITHMS = ["cmaes", "nevergrad", "pybobyqa", "nlopt_direct", "pso", "bayesopt"]


class ObjectiveTrace:
    def __init__(self, problem: ProblemSpec):
        self.problem = problem
        self.nfev = 0
        self.best = np.inf
        self.trace: list[dict[str, float]] = []

    def __call__(self, x: np.ndarray) -> float:
        t0 = time.perf_counter()
        ev = evaluate_theta(x, self.problem)
        self.nfev += 1
        self.best = min(self.best, ev.objective)
        if self.nfev <= 5 or self.nfev % 10 == 0:
            self.trace.append({
                "nfev": float(self.nfev),
                "obj": float(ev.objective),
                "best": float(self.best),
                "feasible": float(ev.feasible),
                "elapsed_eval": float(time.perf_counter() - t0),
            })
        return float(ev.objective)


def available_optional_algorithms() -> dict[str, bool]:
    out = {}
    for name, mod in {
        "cmaes": "cma",
        "nevergrad": "nevergrad",
        "pybobyqa": "pybobyqa",
        "nlopt_direct": "nlopt",
        "pso": "pyswarms",
        "bayesopt": "skopt",
    }.items():
        try:
            __import__(mod)
            out[name] = True
        except Exception:
            out[name] = False
    return out


def scipy_bounds(problem: ProblemSpec, dim: int) -> list[tuple[float, float]]:
    b = problem.parameterization.bounds
    if b is None:
        return [(-0.5, 0.5)] * dim
    return [b] * dim


def optimize_with_algorithm(
    problem: ProblemSpec,
    algorithm_name: str,
    x0: np.ndarray,
    seed: int,
    restart_id: int,
    initialization: str,
    maxiter: int = 300,
    maxfun: int = 2000,
) -> OptimizationResult:
    np.random.seed(seed)
    start = time.perf_counter()
    wrapper = ObjectiveTrace(problem)
    dim = len(np.asarray(x0).reshape(-1))
    bounds = scipy_bounds(problem, dim)
    nit = None
    success = False
    msg = ""

    try:
        if algorithm_name in {"Nelder-Mead", "Powell", "COBYLA", "SLSQP", "L-BFGS-B", "trust-constr"}:
            res = minimize(
                wrapper,
                clip_theta_to_bounds(x0, problem),
                method=algorithm_name,
                bounds=bounds if algorithm_name in {"Powell", "SLSQP", "L-BFGS-B", "trust-constr"} else None,
                options={"maxiter": maxiter},
            )
        elif algorithm_name == "differential_evolution":
            res = differential_evolution(wrapper, bounds=bounds, seed=seed, maxiter=maxiter, polish=False)
        elif algorithm_name == "dual_annealing":
            res = dual_annealing(wrapper, bounds=bounds, seed=seed, maxiter=maxiter)
        elif algorithm_name == "basinhopping":
            minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds}
            res = basinhopping(wrapper, clip_theta_to_bounds(x0, problem), niter=min(100, maxiter), seed=seed, minimizer_kwargs=minimizer_kwargs)
        elif algorithm_name == "shgo":
            res = shgo(wrapper, bounds=bounds, options={"maxev": maxfun})
        else:
            raise ValueError(f"Unsupported algorithm: {algorithm_name}")

        theta = np.asarray(res.x, dtype=float).reshape(-1)
        final = evaluate_theta(theta, problem)
        success = bool(getattr(res, "success", True))
        msg = str(getattr(res, "message", "ok"))
        nit = int(getattr(res, "nit", 0)) if hasattr(res, "nit") else None
        nfev = int(getattr(res, "nfev", wrapper.nfev))
    except Exception as exc:
        theta = np.asarray(x0, dtype=float)
        final = evaluate_theta(theta, problem)
        msg = f"{type(exc).__name__}: {exc}"
        nfev = wrapper.nfev

    runtime = time.perf_counter() - start
    return OptimizationResult(
        algorithm=algorithm_name,
        problem_id=problem.problem_id,
        restart_id=restart_id,
        seed=seed,
        initialization=initialization,
        success=success,
        feasible=bool(final.feasible),
        timeout=False,
        objective=float(final.objective),
        true_objective=float(final.true_objective),
        surrogate_objective=final.surrogate_objective,
        runtime_sec=float(runtime),
        nfev=nfev,
        nit=nit,
        param_norm=float(np.linalg.norm(theta)),
        cond_i_minus_q=final.diagnostics.get("cond_i_minus_q"),
        message=msg,
        theta=theta,
        Q=final.Q,
        diagnostics=final.diagnostics,
        trace=wrapper.trace,
    )
