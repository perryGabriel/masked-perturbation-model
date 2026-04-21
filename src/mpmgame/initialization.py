from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .design_lp import lp_relaxation_design
from .design_projection import projection_design
from .objective_engine import ProblemSpec, q_to_theta, theta_dim, theta_to_q


@dataclass
class InitPoint:
    name: str
    theta: np.ndarray


def _random_theta(problem: ProblemSpec, rng: np.random.Generator, sparse: bool = False) -> np.ndarray:
    d = theta_dim(problem)
    lo, hi = problem.parameterization.bounds or (-0.1, 0.1)
    theta = rng.uniform(lo, hi, size=d)
    if sparse and d > 0:
        keep = rng.random(d) < 0.35
        theta = theta * keep
    return theta


def _dsf_stable_denominator_seed(problem: ProblemSpec) -> np.ndarray | None:
    spec = problem.parameterization
    den_degree = int(spec.q_den_degree)
    den_is_optimized = den_degree > 0 and spec.stability_parameterization != "fixed_den"
    if not den_is_optimized or not spec.shared_denominator:
        return None
    den = np.zeros(den_degree, dtype=float)
    if den_degree > 0:
        # 1 - 0.2 z^-1 (plus optional trailing zero terms) keeps poles strictly inside unit disk.
        den[0] = -0.2
    return den


def _dsf_initializations(problem: ProblemSpec, rng: np.random.Generator, n_random: int) -> list[InitPoint]:
    d = theta_dim(problem)
    spec = problem.parameterization
    num_terms = int(spec.q_num_degree) + 1
    n_free = 0 if num_terms <= 0 else d // num_terms
    den_seed = _dsf_stable_denominator_seed(problem)
    den_count = 0 if den_seed is None else len(den_seed)
    num_count = d - den_count

    inits: list[InitPoint] = []
    base = np.zeros(d, dtype=float)
    if den_seed is not None:
        base[num_count:] = den_seed
    inits.append(InitPoint(name="dsf_near_zero", theta=base))

    for k in range(n_random):
        t = base.copy()
        if num_count > 0 and n_free > 0:
            lo, hi = spec.bounds or (-0.1, 0.1)
            vals = rng.uniform(lo, hi, size=num_count)
            keep = rng.random(n_free) < 0.2
            keep = np.repeat(keep, num_terms)
            t[:num_count] = vals * keep
        inits.append(InitPoint(name=f"dsf_sparse_num_{k}", theta=t))
    return inits


def generate_initializations(problem: ProblemSpec, n_random: int, seed: int) -> list[InitPoint]:
    rng = np.random.default_rng(seed)
    inits: list[InitPoint] = [InitPoint(name="zero", theta=np.zeros(theta_dim(problem)))]

    for k in range(n_random):
        inits.append(InitPoint(name=f"gaussian_{k}", theta=_random_theta(problem, rng, sparse=False)))
        inits.append(InitPoint(name=f"sparse_{k}", theta=_random_theta(problem, rng, sparse=True)))

    if problem.parameterization.kind == "static_hollow":
        # Warm starts from existing baseline solvers
        try:
            pres = projection_design(
                system=problem.system,
                access_model=problem.access_model,
                threat_model=problem.threat_model,
                structure_mask=problem.parameterization.mask,
                max_iter=5,
            )
            inits.append(InitPoint(name="projection_warm", theta=q_to_theta(pres.Q_final, problem)))
        except Exception:
            pass

        try:
            lres = lp_relaxation_design(
                system=problem.system,
                access_model=problem.access_model,
                threat_model=problem.threat_model,
                g_hat=problem.parameterization.g_hat or 0.2,
                structure_mask=problem.parameterization.mask,
            )
            if lres.success:
                inits.append(InitPoint(name="lp_warm", theta=q_to_theta(lres.Q, problem)))
        except Exception:
            pass
    elif problem.parameterization.kind == "dsf_poly":
        try:
            inits.extend(_dsf_initializations(problem, rng=rng, n_random=n_random))
        except Exception:
            inits.append(InitPoint(name="dsf_fallback_zero", theta=np.zeros(theta_dim(problem))))

    return inits


def clip_theta_to_bounds(theta: np.ndarray, problem: ProblemSpec) -> np.ndarray:
    bounds = problem.parameterization.bounds
    if bounds is None:
        return np.asarray(theta, dtype=float)
    lo, hi = bounds
    return np.clip(np.asarray(theta, dtype=float), lo, hi)


def theta_to_q_with_metadata(theta: np.ndarray, problem: ProblemSpec):
    return theta_to_q(theta, problem)
