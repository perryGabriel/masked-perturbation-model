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


def generate_initializations(problem: ProblemSpec, n_random: int, seed: int) -> list[InitPoint]:
    rng = np.random.default_rng(seed)
    inits: list[InitPoint] = [InitPoint(name="zero", theta=np.zeros(theta_dim(problem)))]

    for k in range(n_random):
        inits.append(InitPoint(name=f"gaussian_{k}", theta=_random_theta(problem, rng, sparse=False)))
        inits.append(InitPoint(name=f"sparse_{k}", theta=_random_theta(problem, rng, sparse=True)))

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

    return inits


def clip_theta_to_bounds(theta: np.ndarray, problem: ProblemSpec) -> np.ndarray:
    bounds = problem.parameterization.bounds
    if bounds is None:
        return np.asarray(theta, dtype=float)
    lo, hi = bounds
    return np.clip(np.asarray(theta, dtype=float), lo, hi)


def theta_to_q_with_metadata(theta: np.ndarray, problem: ProblemSpec):
    return theta_to_q(theta, problem)
