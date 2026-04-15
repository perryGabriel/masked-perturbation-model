"""Chapter-2 experiment families and report helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .fmp import ContractSystem, access_matrix, build_contract_system, q_strength_metric, vulnerability_full, vulnerability_single_link
from .bounds import bound_lower_full, bound_lower_single_link, bound_upper_full


@dataclass
class FamilyResult:
    name: str
    data: pd.DataFrame
    notes: list[str]


def _measured(system: ContractSystem, Q: np.ndarray, access_model: str, threat_model: str) -> float:
    pbar = access_matrix(system, Q, access_model)
    if threat_model == "full":
        return vulnerability_full(system, Q, pbar).value
    if threat_model == "single_link":
        return vulnerability_single_link(system, Q, pbar).value
    raise ValueError("threat_model must be 'full' or 'single_link'")


def family_vary_q(system: ContractSystem, Q0: np.ndarray, lambdas: Iterable[float], access_model: str = "w2", threat_model: str = "single_link") -> FamilyResult:
    rows = []
    notes = [f"Threat model: {threat_model}", f"Access model: {access_model}"]
    for lam in lambdas:
        Q = float(lam) * np.asarray(Q0, dtype=float)
        pbar = access_matrix(system, Q, access_model)
        meas = _measured(system, Q, access_model, threat_model)
        b_low_single = bound_lower_single_link(system, Q, pbar).value
        b_low_full = bound_lower_full(system, pbar).value
        b_up_full = bound_upper_full(system, pbar).value
        rows.append({
            "lambda": float(lam),
            "q_norm_2": q_strength_metric(Q),
            "measured": meas,
            "lower_single": b_low_single,
            "lower_full": b_low_full,
            "upper_full": b_up_full,
        })
    return FamilyResult(name="vary_q", data=pd.DataFrame(rows), notes=notes)


def family_vary_alpha(G: np.ndarray, alpha0: np.ndarray, rhos: Iterable[float], Q: np.ndarray, access_model: str = "w2", threat_model: str = "full") -> FamilyResult:
    rows = []
    notes = [f"Threat model: {threat_model}", f"Access model: {access_model}"]
    for rho in rhos:
        alpha = float(rho) * np.asarray(alpha0, dtype=float)
        system = build_contract_system(np.asarray(G, dtype=float), alpha, label=f"rho={rho:.3g}")
        pbar = access_matrix(system, Q, access_model)
        meas = _measured(system, Q, access_model, threat_model)
        lower_full = bound_lower_full(system, pbar).value
        ga_norm = float(np.linalg.svd(system.G @ system.alpha, compute_uv=False)[0])
        rows.append({
            "rho": float(rho),
            "ga_norm": ga_norm,
            "measured": meas,
            "lower_full": lower_full,
        })
    return FamilyResult(name="vary_alpha", data=pd.DataFrame(rows), notes=notes)


def family_vary_access(system: ContractSystem, Q: np.ndarray, access_models: Iterable[str] = ("w1", "w2", "w3"), threat_model: str = "full") -> FamilyResult:
    rows = []
    notes = [f"Threat model: {threat_model}"]
    for am in access_models:
        pbar = access_matrix(system, Q, am)
        meas = _measured(system, Q, am, threat_model)
        lower_full = bound_lower_full(system, pbar).value
        rank_pbar = int(np.linalg.matrix_rank(pbar))
        rows.append({
            "access_model": am,
            "pbar_rank": rank_pbar,
            "measured": meas,
            "lower_full": lower_full,
        })
    return FamilyResult(name="vary_access", data=pd.DataFrame(rows), notes=notes)


def summarize_experiment_results(*families: FamilyResult) -> pd.DataFrame:
    rows = []
    for fam in families:
        d = fam.data
        rows.append({
            "family": fam.name,
            "n_rows": len(d),
            "measured_min": float(d["measured"].min()) if "measured" in d else np.nan,
            "measured_max": float(d["measured"].max()) if "measured" in d else np.nan,
            "notes": " | ".join(fam.notes),
        })
    return pd.DataFrame(rows)


def export_report_markdown(path: str | Path, summary_df: pd.DataFrame, headline: str = "Chapter-2 demo report") -> Path:
    path = Path(path)
    lines = [f"# {headline}", "", "## Summary table", "", summary_df.to_markdown(index=False), ""]
    lines.append("## Interpretation notes")
    lines.append("")
    lines.append("- These are numerical sanity checks on small systems, not theorem proofs.")
    lines.append("- Bounds are compared only within their stated threat/assumption regime.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path
