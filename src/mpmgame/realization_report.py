from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd

from .benchmark_suite import benchmark_problem_registry
from .fmp import access_matrix, make_realization


def _parse_matrix_cell(cell: object) -> np.ndarray:
    if isinstance(cell, np.ndarray):
        return np.asarray(cell, dtype=float)
    return np.asarray(ast.literal_eval(str(cell)), dtype=float)


def load_raw_results(csv_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "Q" in df.columns:
        df["Q"] = df["Q"].map(_parse_matrix_cell)
    if "theta" in df.columns:
        df["theta"] = df["theta"].map(lambda x: np.asarray(ast.literal_eval(str(x)), dtype=float))
    return df


def _problem_by_id(problem_id: str):
    for preset in ("quick", "full"):
        for problem in benchmark_problem_registry(preset):
            if problem.problem_id == problem_id:
                return problem
    raise KeyError(f"Unknown problem_id={problem_id}")


def summarize_row_realization(row: pd.Series) -> dict[str, object]:
    problem = _problem_by_id(str(row["problem_id"]))
    Q = np.asarray(row["Q"], dtype=float)
    P = make_realization(problem.system.G, Q)
    pbar = access_matrix(problem.system, Q, problem.access_model)
    return {
        "problem_id": problem.problem_id,
        "algorithm": str(row["algorithm"]),
        "restart_id": int(row["restart_id"]),
        "objective": float(row["true_objective"]),
        "feasible": bool(row["feasible"]),
        "Q": Q,
        "P": P,
        "Pbar": pbar,
        "I_minus_Q": np.eye(problem.system.n) - Q,
    }


def rows_have_same_contracted_subsystem(row_a: pd.Series, row_b: pd.Series, atol: float = 1e-7) -> bool:
    sa = summarize_row_realization(row_a)
    sb = summarize_row_realization(row_b)
    return np.allclose(sa["P"], sb["P"], atol=atol, rtol=0.0)


def realization_markdown(summary: dict[str, object], include_sympy: bool = True) -> str:
    lines = [
        f"# Realization inspection: {summary['problem_id']} / {summary['algorithm']} / restart {summary['restart_id']}",
        "",
        f"- Feasible: `{summary['feasible']}`",
        f"- Objective: `{summary['objective']:.6g}`",
        "",
        "## Numerical matrices",
        "",
        "### Q",
        "```text",
        np.array2string(summary["Q"], precision=5, suppress_small=True),
        "```",
        "",
        "### Contracted subsystem P=(I-Q)G",
        "```text",
        np.array2string(summary["P"], precision=5, suppress_small=True),
        "```",
    ]
    if include_sympy:
        try:
            import sympy as sp

            q_sp = sp.Matrix(np.asarray(summary["Q"], dtype=float))
            p_sp = sp.Matrix(np.asarray(summary["P"], dtype=float))
            lines += [
                "",
                "## SymPy forms",
                "",
                "```text",
                f"Q_sympy = {sp.srepr(q_sp)}",
                f"P_sympy = {sp.srepr(p_sp)}",
                "```",
            ]
        except Exception:
            lines += ["", "_SymPy not available in this environment; numeric form shown above._"]
    return "\n".join(lines) + "\n"


def write_realization_markdown_from_csv(
    csv_path: str | Path,
    out_path: str | Path,
    problem_id: str,
    algorithm: str | None = None,
) -> Path:
    df = load_raw_results(csv_path)
    scope = df[(df["problem_id"] == problem_id) & (df["feasible"] == True)]  # noqa: E712
    if algorithm is not None:
        scope = scope[scope["algorithm"] == algorithm]
    if scope.empty:
        raise ValueError("No feasible rows available for the requested filters.")
    row = scope.loc[scope["true_objective"].idxmin()]
    text = realization_markdown(summarize_row_realization(row))
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return out
