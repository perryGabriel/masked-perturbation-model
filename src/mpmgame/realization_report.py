from __future__ import annotations

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .benchmark_suite import benchmark_problem_registry
from .dsf_param import DynamicQ, deserialize_dynamic_q_from_row
from .fmp import access_matrix, make_realization_checked


def _parse_matrix_cell(cell: object) -> np.ndarray:
    if isinstance(cell, np.ndarray):
        return np.asarray(cell, dtype=float)
    return np.asarray(ast.literal_eval(str(cell)), dtype=float)


def _parse_optional_json_cell(cell: object) -> object | None:
    if cell is None:
        return None
    if isinstance(cell, (dict, list)):
        return cell
    text = str(cell).strip()
    if text in {"", "nan", "None"}:
        return None
    try:
        return json.loads(text)
    except Exception:
        try:
            return ast.literal_eval(text)
        except Exception:
            return None


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


def _dynamic_q_from_row(row: pd.Series) -> DynamicQ | None:
    payload_col = row.get("dynamic_q_payload", None)
    if payload_col is not None and str(payload_col).strip() not in {"", "nan", "None"}:
        payload = payload_col if isinstance(payload_col, dict) else json.loads(str(payload_col))
        return deserialize_dynamic_q_from_row(payload)

    q_coeffs = _parse_optional_json_cell(row.get("Q_coeffs", None))
    if isinstance(q_coeffs, dict) and q_coeffs.get("num_coeffs", None) is not None:
        payload = {
            "kind": "dsf_dynamic_q",
            "num_coeffs": q_coeffs["num_coeffs"],
            "den_coeffs": q_coeffs.get("den_coeffs", None),
            "metadata": _parse_optional_json_cell(row.get("Q_meta", None)) or {},
        }
        return deserialize_dynamic_q_from_row(payload)
    return None


def _stable_denominator_flags(dynamic_q: DynamicQ) -> dict[str, object]:
    den = dynamic_q.den_coeffs
    if den is None or np.asarray(den).size == 0:
        return {"stable": True, "max_root_magnitude": 0.0, "unstable_links": 0}

    if den.ndim == 1:
        roots = np.roots(np.concatenate([[1.0], den]))
        max_mag = float(np.max(np.abs(roots))) if roots.size else 0.0
        return {"stable": bool(max_mag < 1.0), "max_root_magnitude": max_mag, "unstable_links": int(max_mag >= 1.0)}

    unstable = 0
    max_mag = 0.0
    n = den.shape[0]
    for i in range(n):
        for j in range(n):
            roots = np.roots(np.concatenate([[1.0], den[i, j, :]]))
            entry_max = float(np.max(np.abs(roots))) if roots.size else 0.0
            max_mag = max(max_mag, entry_max)
            if entry_max >= 1.0:
                unstable += 1
    return {"stable": unstable == 0, "max_root_magnitude": max_mag, "unstable_links": unstable}


def summarize_row_realization(row: pd.Series) -> dict[str, object]:
    problem = _problem_by_id(str(row["problem_id"]))
    Q = np.asarray(row["Q"], dtype=float)
    param_kind = str(row.get("param_kind", "static_hollow") or "static_hollow")

    q_for_realization: np.ndarray = Q
    dsf_summary: dict[str, object] | None = None
    symbolic_forms: dict[str, str] | None = None

    dynamic_q = _dynamic_q_from_row(row)
    if dynamic_q is not None:
        param_kind = "dsf_poly"
        q_for_realization = dynamic_q.to_tf_matrix()

        coeff_abs = np.abs(dynamic_q.num_coeffs)
        nonzero_links = int(np.count_nonzero(np.any(coeff_abs > 1e-12, axis=2)))
        den_flags = _stable_denominator_flags(dynamic_q)
        dsf_summary = {
            "num_degree": int(dynamic_q.num_degree),
            "den_degree": int(dynamic_q.den_degree),
            "nonzero_links": nonzero_links,
            "stable_denominator": bool(den_flags["stable"]),
            "max_root_magnitude": float(den_flags["max_root_magnitude"]),
            "unstable_links": int(den_flags["unstable_links"]),
        }

        q_tf = dynamic_q.to_tf_matrix()
        n = q_tf.shape[0]
        rendered: list[list[str]] = []
        for i in range(n):
            row_forms: list[str] = []
            for j in range(n):
                tf_ij = q_tf[i, j]
                num = np.asarray(tf_ij.num[0][0], dtype=float)
                den = np.asarray(tf_ij.den[0][0], dtype=float)
                row_forms.append(f"({np.array2string(num, precision=4)})/({np.array2string(den, precision=4)})")
            rendered.append(row_forms)
        symbolic_forms = {"Q_tf": str(rendered)}

    try:
        P = make_realization_checked(problem.system.G, q_for_realization)
        pbar = access_matrix(problem.system, q_for_realization, problem.access_model)
    except Exception:
        P = make_realization_checked(problem.system.G, Q)
        pbar = access_matrix(problem.system, Q, problem.access_model)
    return {
        "problem_id": problem.problem_id,
        "algorithm": str(row["algorithm"]),
        "restart_id": int(row["restart_id"]),
        "objective": float(row["true_objective"]),
        "feasible": bool(row["feasible"]),
        "param_kind": param_kind,
        "Q": Q,
        "P": P,
        "Pbar": pbar,
        "I_minus_Q": np.eye(problem.system.n) - Q,
        "dsf_summary": dsf_summary,
        "symbolic_forms": symbolic_forms,
    }


def rows_have_same_contracted_subsystem(row_a: pd.Series, row_b: pd.Series, atol: float = 1e-7) -> bool:
    sa = summarize_row_realization(row_a)
    sb = summarize_row_realization(row_b)
    return np.allclose(sa["P"], sb["P"], atol=atol, rtol=0.0)


def realization_markdown(summary: dict[str, object], include_sympy: bool = True) -> str:
    lines = [
        f"# Realization inspection: {summary['problem_id']} / {summary['algorithm']} / restart {summary['restart_id']}",
        "",
        f"- Parameterization: `{summary.get('param_kind', 'static_hollow')}`",
        f"- Feasible: `{summary['feasible']}`",
        f"- Objective: `{summary['objective']:.6g}`",
    ]

    dsf_summary = summary.get("dsf_summary", None)
    if isinstance(dsf_summary, dict):
        lines += [
            "",
            "## DSF summary",
            "",
            f"- Degrees: numerator={dsf_summary['num_degree']}, denominator={dsf_summary['den_degree']}",
            f"- Nonzero links: {dsf_summary['nonzero_links']}",
            (
                f"- Stability: stable_denominator={dsf_summary['stable_denominator']}, "
                f"max_root_magnitude={dsf_summary['max_root_magnitude']:.6g}, "
                f"unstable_links={dsf_summary['unstable_links']}"
            ),
        ]
        symbolic_forms = summary.get("symbolic_forms", None)
        if isinstance(symbolic_forms, dict) and symbolic_forms.get("Q_tf", None):
            lines += ["", "### Symbolic / transfer-form entries", "```text", str(symbolic_forms["Q_tf"]), "```"]

    lines += [
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
