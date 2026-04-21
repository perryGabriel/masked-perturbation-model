from __future__ import annotations

import numpy as np

from mpmgame.benchmark_suite import benchmark_problem_registry
from mpmgame.fmp import build_contract_system, make_realization_checked
from mpmgame.objective_engine import (
    ParameterizationSpec,
    ProblemSpec,
    evaluate_theta,
    q_to_theta,
    theta_dim,
    theta_to_q,
)
from mpmgame.tf_tools import tf_matrix


def _base_problem() -> ProblemSpec:
    return benchmark_problem_registry("quick")[0]


def test_static_hollow_regression_small_fixed_cases() -> None:
    """Regression anchor: static_hollow outputs for fixed tiny benchmark cases."""
    expected = {
        "toy2_w1_full": 0.7024908976438582,
        "toy3_w2_single": 1.0052924849143805,
        "toy3_w3_full": 0.07689330141592574,
    }
    for problem in benchmark_problem_registry("quick"):
        d = theta_dim(problem)
        theta = np.linspace(-0.05, 0.05, d)
        ev = evaluate_theta(theta, problem)
        assert ev.feasible
        np.testing.assert_allclose(ev.true_objective, expected[problem.problem_id], rtol=0.0, atol=1e-12)


def test_dsf_poly_theta_dimension_and_mapping_roundtrip() -> None:
    base = _base_problem()
    dsf_problem = ProblemSpec(
        problem_id="dsf_mapping",
        system=base.system,
        access_model=base.access_model,
        threat_model=base.threat_model,
        parameterization=ParameterizationSpec(kind="dsf_poly", q_num_degree=0, q_den_degree=0, zero_diagonal=True),
    )

    assert theta_dim(dsf_problem) == 2
    theta = np.array([0.07, -0.03])
    q_static, _ = theta_to_q(theta, dsf_problem)
    back = q_to_theta(q_static, dsf_problem)
    np.testing.assert_allclose(theta, back)


def test_dsf_poly_zero_diagonal_enforcement() -> None:
    base = _base_problem()
    dsf_problem = ProblemSpec(
        problem_id="dsf_zero_diag",
        system=base.system,
        access_model=base.access_model,
        threat_model=base.threat_model,
        parameterization=ParameterizationSpec(kind="dsf_poly", q_num_degree=1, q_den_degree=1, zero_diagonal=True),
    )
    theta = np.array([0.1, 0.2, -0.3, 0.4])
    q_static, meta = theta_to_q(theta, dsf_problem)

    assert np.allclose(np.diag(q_static), 0.0)
    num_coeffs = np.asarray(meta["num_coeffs"])
    assert np.allclose(num_coeffs[0, 0, :], 0.0)
    assert np.allclose(num_coeffs[1, 1, :], 0.0)


def test_dsf_poly_properness_and_stability_flags() -> None:
    base = _base_problem()

    improper_problem = ProblemSpec(
        problem_id="dsf_improper",
        system=base.system,
        access_model=base.access_model,
        threat_model=base.threat_model,
        parameterization=ParameterizationSpec(
            kind="dsf_poly", q_num_degree=1, q_den_degree=0, stability_parameterization="fixed_den", zero_diagonal=True
        ),
    )
    improper = evaluate_theta(np.array([0.2, 0.3, -0.1, 0.4]), improper_problem)
    assert not improper.feasible
    assert "improper_entry" in improper.diagnostics["reasons"]

    unstable_problem = ProblemSpec(
        problem_id="dsf_unstable",
        system=base.system,
        access_model=base.access_model,
        threat_model=base.threat_model,
        parameterization=ParameterizationSpec(
            kind="dsf_poly",
            q_num_degree=1,
            q_den_degree=1,
            stability_parameterization="free",
            shared_denominator=True,
            zero_diagonal=True,
        ),
    )
    unstable = evaluate_theta(np.array([0.05, 0.01, -0.02, 0.02, -1.2]), unstable_problem)
    assert not unstable.feasible
    assert "unstable_denominator" in unstable.diagnostics["reasons"]


def test_dsf_poly_dynamic_p_construction_sanity() -> None:
    G = np.diag([0.5, 0.8])
    alpha = np.array([[0.0, 0.05], [0.03, 0.0]])
    system = build_contract_system(G, alpha, label="dyn_p")
    _ = system  # explicit tiny fixture creation for readability

    num = np.zeros((2, 2, 2), dtype=float)
    num[0, 1, :] = [0.2, 0.1]
    num[1, 0, :] = [-0.05, 0.02]
    den = np.array([0.3], dtype=float)

    from mpmgame.dsf_param import DynamicQ

    dynamic_q = DynamicQ.from_tensors(num_coeffs=num, den_coeffs=den, zero_diagonal=True)
    p_tf = make_realization_checked(tf_matrix(G.tolist()), dynamic_q.to_tf_matrix())

    assert p_tf.shape == (2, 2)
    assert len(p_tf[0, 1].num[0][0]) >= 2
    assert len(p_tf[0, 1].den[0][0]) >= 2


def test_dsf_poly_objective_tiny_frequency_grid() -> None:
    base = _base_problem()
    freq_grid = np.array([0.0, 0.5, 1.0], dtype=float)
    dsf_problem = ProblemSpec(
        problem_id="dsf_freq_obj",
        system=base.system,
        access_model=base.access_model,
        threat_model=base.threat_model,
        parameterization=ParameterizationSpec(kind="dsf_poly", q_num_degree=0, q_den_degree=0, zero_diagonal=True),
        freq_grid=freq_grid,
    )

    ev = evaluate_theta(np.array([0.05, -0.03]), dsf_problem)
    assert ev.feasible
    assert np.isfinite(ev.objective)
    assert len(ev.diagnostics["cond_i_minus_q_grid"]) == len(freq_grid)
