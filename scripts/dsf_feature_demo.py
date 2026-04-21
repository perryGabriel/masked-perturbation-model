#!/usr/bin/env python3
"""End-to-end demo of new DSF functionality.

This script is intentionally explicit about assumptions so results are easy to
interpret when compared against canonical DSF definitions (e.g., Warnick,
Goncalves, et al.).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mpmgame.benchmark_suite import benchmark_problem_registry, run_benchmark_suite
from mpmgame.dsf_param import DynamicQ
from mpmgame.fmp import make_realization_dynamic, vulnerability_full, vulnerability_single_link
from mpmgame.objective_engine import ParameterizationSpec, ProblemSpec, evaluate_theta, theta_dim, theta_to_q
from mpmgame.realization_report import write_realization_markdown_from_csv
from mpmgame.tf_tools import tf_matrix


def _section(title: str) -> None:
    print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")


def _demo_static_compatibility() -> None:
    _section("1) Legacy static path remains unchanged")
    static_problem = benchmark_problem_registry("quick")[0]
    d = theta_dim(static_problem)
    theta0 = np.zeros(d)
    ev = evaluate_theta(theta0, static_problem)
    print(f"problem_id={static_problem.problem_id}")
    print(f"kind={static_problem.parameterization.kind}, theta_dim={d}")
    print(f"feasible={ev.feasible}, objective={ev.objective:.6g}, true_objective={ev.true_objective:.6g}")


def _build_demo_dsf_problem() -> ProblemSpec:
    static_problem = benchmark_problem_registry("quick")[0]
    param = ParameterizationSpec(
        kind="dsf_poly",
        bounds=(-0.2, 0.2),
        g_hat=0.6,
        q_num_degree=1,
        q_den_degree=1,
        shared_denominator=True,
        stability_parameterization="fixed_den",
        zero_diagonal=True,
    )
    return ProblemSpec(
        problem_id="demo_dsf_qp",
        system=static_problem.system,
        access_model="w1",
        threat_model="full",
        parameterization=param,
        freq_grid=np.logspace(-2, 2, 120),
        cond_max=1e6,
    )


def _demo_dynamic_qp(problem: ProblemSpec) -> None:
    _section("2) DSF parameterization, Q(z) construction, and P=(I-Q)G")
    d = theta_dim(problem)
    theta = np.linspace(-0.05, 0.05, d)
    q_static, meta = theta_to_q(theta, problem)

    dynamic_q = DynamicQ.from_tensors(
        num_coeffs=np.asarray(meta["num_coeffs"], dtype=float),
        den_coeffs=meta.get("den_coeffs", None),
        mask=problem.parameterization.mask,
        zero_diagonal=problem.parameterization.zero_diagonal,
        metadata={"source": "dsf_feature_demo"},
    )

    violations = dynamic_q.verify_constraints()
    q_tf = dynamic_q.to_tf_matrix()
    g_tf = tf_matrix(problem.system.G.tolist())
    p_tf = make_realization_dynamic(g_tf, q_tf)

    print(f"problem_id={problem.problem_id}")
    print(f"theta_dim={d}, q_num_degree={problem.parameterization.q_num_degree}, q_den_degree={problem.parameterization.q_den_degree}")
    print(f"Q(0) offdiag snapshot:\n{q_static}")
    print(f"constraint_violations={violations if violations else 'none'}")
    print(f"dynamic P shape={p_tf.shape}, Q shape={q_tf.shape}")

    # For the dynamic demo we use w1 semantics directly: Pbar = (I-Q)G.
    pbar = p_tf
    vf = vulnerability_full(problem.system, q_tf, pbar, freq_grid=problem.freq_grid)
    vs = vulnerability_single_link(problem.system, q_tf, pbar, freq_grid=problem.freq_grid)

    print("dynamic vulnerability estimates (frequency-grid approximations):")
    print(f"  full attacker: value={vf.value:.6g}, approximate={vf.approximate}, diagnostics={vf.diagnostics}")
    print(f"  single-link : value={vs.value:.6g}, approximate={vs.approximate}, diagnostics={vs.diagnostics}")


def _demo_serialization_and_report(outdir: Path) -> None:
    _section("3) Mixed benchmark run + DSF serialization/report")
    static_problem = benchmark_problem_registry("quick")[0]

    dsf_problem = ProblemSpec(
        problem_id="toy2_w1_full_dsf",
        system=static_problem.system,
        access_model=static_problem.access_model,
        threat_model=static_problem.threat_model,
        parameterization=ParameterizationSpec(
            kind="dsf_poly",
            bounds=(-0.2, 0.2),
            g_hat=0.6,
            q_num_degree=0,
            q_den_degree=0,
            zero_diagonal=True,
        ),
    )

    results_dir = outdir / "results"
    reports_dir = outdir / "reports"
    raw_df, summary_df = run_benchmark_suite(
        problem_specs=[static_problem, dsf_problem],
        algorithms=["Powell"],
        n_restarts=1,
        timeout_per_run_sec=10,
        output_dir=results_dir,
        report_dir=reports_dir,
        seed=13,
        show_progress=False,
        maxiter=20,
    )

    print(f"raw_rows={len(raw_df)}, summary_rows={len(summary_df)}")
    for col in ["problem_id", "algorithm", "feasible", "param_kind", "Q_coeffs", "Q_meta"]:
        print(f"column_present[{col}]={col in raw_df.columns}")

    dsf_rows = raw_df[raw_df["problem_id"] == "toy2_w1_full_dsf"]
    any_payload = bool(dsf_rows["Q_coeffs"].notna().any()) if not dsf_rows.empty else False
    print(f"dsf_rows={len(dsf_rows)}, any_dsf_payload={any_payload}")

    out_md = reports_dir / "demo_static_realization.md"
    write_realization_markdown_from_csv(
        results_dir / "benchmark_raw_results.csv",
        out_md,
        problem_id="toy2_w1_full",
    )
    print(f"wrote realization report (static compatibility check): {out_md}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DSF feature demonstration.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/dsf_feature_demo"),
        help="Output directory for generated benchmark/report artifacts.",
    )
    parser.add_argument(
        "--skip-benchmark",
        action="store_true",
        help="Skip the benchmark/report phase and only run static+dynamic API demos.",
    )
    args = parser.parse_args()

    _demo_static_compatibility()
    dsf_problem = _build_demo_dsf_problem()
    _demo_dynamic_qp(dsf_problem)

    if not args.skip_benchmark:
        _demo_serialization_and_report(args.outdir)


if __name__ == "__main__":
    main()
