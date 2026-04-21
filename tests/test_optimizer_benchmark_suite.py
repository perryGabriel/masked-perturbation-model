import time

import numpy as np

from mpmgame.benchmark_suite import benchmark_problem_registry, run_benchmark_suite
from mpmgame.objective_engine import ParameterizationSpec, ProblemSpec, evaluate_theta, q_to_theta, theta_dim, theta_to_q
from mpmgame.realization_report import (
    load_raw_results,
    rows_have_same_contracted_subsystem,
    summarize_row_realization,
    write_realization_markdown_from_csv,
)
from mpmgame.timeouts import run_with_timeout


def _sleepy(x):
    time.sleep(0.3)
    return x


def _small_problem():
    return benchmark_problem_registry("quick")[0]


def test_timeout_enforcement():
    out = run_with_timeout(_sleepy, 0.05, 1)
    assert out.timed_out


def test_parameter_vector_mapping_roundtrip():
    problem = _small_problem()
    d = theta_dim(problem)
    theta = np.linspace(-0.05, 0.05, d)
    Q, _ = theta_to_q(theta, problem)
    back = q_to_theta(Q, problem)
    np.testing.assert_allclose(theta, back)


def test_feasibility_check_catches_bound_violation():
    problem = _small_problem()
    d = theta_dim(problem)
    theta = np.zeros(d)
    assert evaluate_theta(theta, problem).feasible

    bad_problem = ProblemSpec(
        problem_id=problem.problem_id,
        system=problem.system,
        access_model=problem.access_model,
        threat_model=problem.threat_model,
        parameterization=ParameterizationSpec(kind="static_hollow", bounds=(-0.01, 0.01), g_hat=0.02),
    )
    bad_theta = np.full(d, 0.1)
    bad = evaluate_theta(bad_theta, bad_problem)
    assert not bad.feasible


def test_aggregation_and_seed_reproducibility(tmp_path):
    probs = benchmark_problem_registry("quick")[:1]
    kwargs = dict(
        problem_specs=probs,
        algorithms=["Nelder-Mead"],
        n_restarts=2,
        timeout_per_run_sec=5,
        output_dir=tmp_path / "r1",
        report_dir=tmp_path / "p1",
        seed=123,
        show_progress=False,
        maxiter=20,
    )
    raw1, sum1 = run_benchmark_suite(**kwargs)
    raw2, sum2 = run_benchmark_suite(**{**kwargs, "output_dir": tmp_path / "r2", "report_dir": tmp_path / "p2"})
    assert not sum1.empty
    a1 = raw1[(raw1["algorithm"] == "Nelder-Mead") & (raw1["restart_id"] == 0)]["true_objective"].iloc[0]
    a2 = raw2[(raw2["algorithm"] == "Nelder-Mead") & (raw2["restart_id"] == 0)]["true_objective"].iloc[0]
    assert np.isclose(a1, a2)


def test_optional_optimizer_unavailable_is_graceful(tmp_path):
    raw, _ = run_benchmark_suite(
        problem_specs=benchmark_problem_registry("quick")[:1],
        algorithms=["nevergrad"],
        n_restarts=1,
        timeout_per_run_sec=2,
        output_dir=tmp_path / "rr",
        report_dir=tmp_path / "pp",
        show_progress=False,
    )
    never = raw[raw["algorithm"] == "nevergrad"]
    assert not never.empty
    assert never["message"].str.contains("unavailable").any()


def test_tqdm_non_interactive(tmp_path):
    raw, _ = run_benchmark_suite(
        problem_specs=benchmark_problem_registry("quick")[:1],
        algorithms=["Powell"],
        n_restarts=1,
        timeout_per_run_sec=3,
        output_dir=tmp_path / "r",
        report_dir=tmp_path / "p",
        show_progress=False,
        maxiter=10,
    )
    assert not raw.empty


def test_incremental_notes_and_realization_export(tmp_path):
    out_dir = tmp_path / "r"
    rep_dir = tmp_path / "p"
    notes_dir = tmp_path / "notes"
    raw, _ = run_benchmark_suite(
        problem_specs=benchmark_problem_registry("quick")[:1],
        algorithms=["Powell"],
        n_restarts=2,
        timeout_per_run_sec=3,
        output_dir=out_dir,
        report_dir=rep_dir,
        incremental_notes_dir=notes_dir,
        show_progress=False,
        maxiter=10,
    )
    assert not raw.empty
    note_files = list(notes_dir.rglob("*.md"))
    assert note_files
    text = note_files[0].read_text()
    assert "Failure modes observed" in text

    csv_path = out_dir / "benchmark_raw_results.csv"
    loaded = load_raw_results(csv_path)
    feasible = loaded[(loaded["problem_id"] == "toy2_w1_full") & (loaded["feasible"] == True)]  # noqa: E712
    assert not feasible.empty
    row = feasible.iloc[0]
    summary = summarize_row_realization(row)
    assert summary["Q"].shape == (2, 2)

    md_out = tmp_path / "realization.md"
    write_realization_markdown_from_csv(csv_path, md_out, problem_id="toy2_w1_full")
    assert md_out.exists()
    assert "Contracted subsystem" in md_out.read_text()

    assert rows_have_same_contracted_subsystem(row, row)
