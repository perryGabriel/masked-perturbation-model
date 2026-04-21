from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def aggregate_results(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    g = df.groupby(["problem_id", "algorithm"], dropna=False)
    rows = []
    for (pid, algo), d in g:
        feas = d[d["feasible"] == True]  # noqa: E712
        best_val = float(feas["true_objective"].min()) if not feas.empty else np.inf
        median_val = float(feas["true_objective"].median()) if not feas.empty else np.inf
        best_runtime = float(feas["runtime_sec"].min()) if not feas.empty else np.inf
        frac_feasible = float((d["feasible"] == True).mean())  # noqa: E712
        frac_timeout = float((d["timeout"] == True).mean())  # noqa: E712
        robustness = 0.7 * frac_feasible + 0.3 * (1.0 / (1.0 + best_val if np.isfinite(best_val) else 1e6))
        rows.append(
            {
                "problem_id": pid,
                "algorithm": algo,
                "runs": len(d),
                "best_feasible_objective": best_val,
                "median_feasible_objective": median_val,
                "fraction_feasible": frac_feasible,
                "fraction_timeout": frac_timeout,
                "mean_runtime_sec": float(d["runtime_sec"].mean()),
                "best_runtime_to_feasible": best_runtime,
                "robustness_score": float(robustness),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["problem_id", "best_feasible_objective", "robustness_score"], ascending=[True, True, False])


def write_report_markdown(report_path: Path, raw_df: pd.DataFrame, summary_df: pd.DataFrame, setup_text: str = "") -> Path:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Nonlinear Optimizer Benchmark Report", "", "## Setup", "", setup_text, ""]
    try:
        sample_table = raw_df.head(30).to_markdown(index=False)
        summary_table = summary_df.to_markdown(index=False)
    except Exception:
        sample_table = raw_df.head(30).to_string(index=False)
        summary_table = summary_df.to_string(index=False)
    lines += ["## Per-run sample", "", sample_table, ""]
    lines += ["## Aggregated summary", "", summary_table, ""]
    if not summary_df.empty:
        best = summary_df.loc[summary_df.groupby("problem_id")["best_feasible_objective"].idxmin()][["problem_id", "algorithm", "best_feasible_objective"]]
        try:
            best_table = best.to_markdown(index=False)
        except Exception:
            best_table = best.to_string(index=False)
        lines += ["## Best algorithm per problem", "", best_table, ""]
    report_path.write_text("\n".join(lines))
    return report_path


def _failure_summary(run_df: pd.DataFrame) -> str:
    failures = run_df[run_df["feasible"] != True]  # noqa: E712
    if failures.empty:
        return "No infeasible runs for this algorithm/problem slice."
    counts = failures["message"].fillna("unknown").value_counts().head(5)
    return "; ".join([f"{msg}: {int(ct)}" for msg, ct in counts.items()])


def write_incremental_result_markdown(
    report_path: Path,
    run_df: pd.DataFrame,
    problem_id: str,
    algorithm: str,
) -> Path:
    """Write per-problem/per-algorithm progress notes.

    This file is intended to be updated after each restart so interpretation can
    keep pace with data generation.
    """
    report_path.parent.mkdir(parents=True, exist_ok=True)
    scope = run_df[(run_df["problem_id"] == problem_id) & (run_df["algorithm"] == algorithm)].copy()
    lines = [
        f"# Optimization Notes: {problem_id} / {algorithm}",
        "",
        f"- Total runs so far: {len(scope)}",
    ]
    feasible = scope[scope["feasible"] == True]  # noqa: E712
    if feasible.empty:
        lines.append("- Feasible runs so far: 0")
        lines.append("- Best feasible objective: n/a")
    else:
        best_idx = feasible["true_objective"].idxmin()
        best = feasible.loc[best_idx]
        lines.append(f"- Feasible runs so far: {len(feasible)}")
        lines.append(f"- Best feasible objective: {best['true_objective']:.6g}")
        lines.append(f"- Best restart id: {int(best['restart_id'])}")
        lines.append(f"- Best runtime (sec): {float(best['runtime_sec']):.4g}")
    lines += [
        "",
        "## Failure modes observed",
        "",
        _failure_summary(scope),
        "",
        "## Recent runs",
        "",
    ]
    try:
        table = scope.sort_values("restart_id").tail(10)[
            ["restart_id", "seed", "initialization", "success", "feasible", "timeout", "true_objective", "runtime_sec", "message"]
        ].to_markdown(index=False)
    except Exception:
        table = scope.sort_values("restart_id").tail(10).to_string(index=False)
    lines.append(table)
    report_path.write_text("\n".join(lines))
    return report_path


def generate_plots(raw_df: pd.DataFrame, summary_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    if raw_df.empty:
        return paths

    # 1) best feasible objective by algorithm
    b = summary_df.groupby("algorithm")["best_feasible_objective"].min().sort_values()
    fig, ax = plt.subplots(figsize=(8, 4))
    b.plot(kind="bar", ax=ax)
    ax.set_ylabel("best feasible objective")
    p = out_dir / "best_feasible_bar.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 2) boxplot of objective over restarts
    fig, ax = plt.subplots(figsize=(9, 4))
    raw_df.boxplot(column="true_objective", by="algorithm", ax=ax, rot=45)
    ax.set_title("Objective across restarts")
    fig.suptitle("")
    p = out_dir / "objective_boxplot.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 3) runtime vs best-value scatter
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(raw_df["runtime_sec"], raw_df["true_objective"], alpha=0.7)
    ax.set_xlabel("runtime (s)"); ax.set_ylabel("true objective")
    p = out_dir / "runtime_vs_objective.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 4) feasibility rates
    f = raw_df.groupby("algorithm")["feasible"].mean().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(8, 4)); f.plot(kind="bar", ax=ax)
    ax.set_ylabel("feasibility rate")
    p = out_dir / "feasibility_rate.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 5) timeout counts
    t = raw_df.groupby("algorithm")["timeout"].sum().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(8, 4)); t.plot(kind="bar", ax=ax)
    ax.set_ylabel("timeout count")
    p = out_dir / "timeout_count.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 6) convergence traces (first feasible per algorithm)
    fig, ax = plt.subplots(figsize=(8, 4))
    for algo, d in raw_df.groupby("algorithm"):
        trace = d.iloc[0].get("trace")
        if isinstance(trace, list) and trace:
            x = [z["nfev"] for z in trace]
            y = [z["best"] for z in trace]
            ax.plot(x, y, label=algo)
    if ax.lines:
        ax.legend(fontsize=8)
    ax.set_xlabel("nfev"); ax.set_ylabel("best-so-far objective")
    p = out_dir / "convergence_traces.png"
    fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 7) best Q heatmap
    feas = raw_df[raw_df["feasible"] == True]  # noqa: E712
    if not feas.empty:
        best_row = feas.loc[feas["true_objective"].idxmin()]
        Q = np.array(best_row["Q"])
        fig, ax = plt.subplots(figsize=(4, 4))
        im = ax.imshow(Q, cmap="coolwarm")
        fig.colorbar(im, ax=ax)
        ax.set_title(f"Best Q ({best_row['algorithm']})")
        p = out_dir / "best_q_heatmap.png"
        fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    # 8) baseline comparison
    baselines = raw_df[raw_df["algorithm"].isin(["baseline_zero", "baseline_projection", "baseline_lp"])]
    if not baselines.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        baselines.groupby("algorithm")["true_objective"].min().plot(kind="bar", ax=ax)
        ax.set_ylabel("best objective")
        p = out_dir / "baseline_comparison.png"
        fig.tight_layout(); fig.savefig(p); plt.close(fig); paths.append(p)

    return paths
