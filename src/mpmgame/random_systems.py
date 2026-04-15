"""Random small-system generators and bound experiment helpers."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .fmp import build_contract_system, access_matrix, q_strength_metric, vulnerability_full, vulnerability_single_link, well_posed
from .bounds import bound_lower_full, bound_lower_single_link, bound_upper_full


@dataclass
class RandomExperimentResult:
    accepted: pd.DataFrame
    rejected: pd.DataFrame


def random_static_system(n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate conservative small static system (diagonal-dominant G, hollow alpha, hollow Q)."""
    g_diag = rng.uniform(0.1, 0.9, size=n)
    G = np.diag(g_diag)
    alpha = rng.uniform(-0.4, 0.4, size=(n, n))
    np.fill_diagonal(alpha, 0.0)
    Q = rng.uniform(-0.2, 0.2, size=(n, n))
    np.fill_diagonal(Q, 0.0)
    return G, alpha, Q


def random_bound_experiments(
    n_samples: int = 100,
    n: int = 3,
    seed: int = 0,
    access_model: str = "w2",
    threat_model: str = "single_link",
    cond_max: float = 1e6,
) -> RandomExperimentResult:
    rng = np.random.default_rng(seed)
    acc = []
    rej = []
    for k in range(n_samples):
        G, alpha, Q = random_static_system(n=n, rng=rng)
        system = build_contract_system(G, alpha, label=f"rand_{k}")
        if not well_posed(system, Q, cond_max=cond_max):
            rej.append({"sample": k, "reason": "well_posedness"})
            continue
        pbar = access_matrix(system, Q, access_model)

        if threat_model == "single_link":
            measured = vulnerability_single_link(system, Q, pbar).value
            lower = bound_lower_single_link(system, Q, pbar).value
            upper = np.nan
        elif threat_model == "full":
            measured = vulnerability_full(system, Q, pbar).value
            lower = bound_lower_full(system, pbar).value
            upper = bound_upper_full(system, pbar).value
        else:
            raise ValueError("threat_model must be 'single_link' or 'full'")

        ratio_low = measured / lower if lower > 0 else np.nan
        ratio_up = measured / upper if (not np.isnan(upper) and upper > 0) else np.nan
        acc.append({
            "sample": k,
            "measured": measured,
            "lower_bound": lower,
            "upper_bound": upper,
            "ratio_to_lower": ratio_low,
            "ratio_to_upper": ratio_up,
            "q_norm": q_strength_metric(Q),
        })

    return RandomExperimentResult(accepted=pd.DataFrame(acc), rejected=pd.DataFrame(rej))
