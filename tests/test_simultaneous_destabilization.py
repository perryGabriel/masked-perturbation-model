import numpy as np

import mpmgame as mpm


def _toy_problem(seed: int = 17) -> mpm.SimultaneousProblem:
    rng = np.random.default_rng(seed)
    masks = [
        mpm.MaskDefense("∇1", np.array([1, 1, 0, 1, 0, 0], dtype=int), threshold=1.9),
        mpm.MaskDefense("∇2", np.array([1, 0, 1, 0, 1, 0], dtype=int), threshold=1.7),
        mpm.MaskDefense("∇3", np.array([0, 1, 1, 0, 0, 1], dtype=int), threshold=1.6),
    ]
    attacks = [
        mpm.AttackCandidate(f"Δ{i+1}", rng.uniform(0.2, 1.2, size=6))
        for i in range(3)
    ]
    sensitivity = np.array([1.0, 1.2, 0.8, 1.1, 1.05, 0.95], dtype=float)
    return mpm.SimultaneousProblem(masks=masks, attacks=attacks, sensitivity=sensitivity)


def test_matrix_shapes_are_correct():
    problem = _toy_problem()
    result = mpm.analyze_simultaneous_attacks(problem, include_additive=True, additive_max_terms=2)

    assert result.cross_success.shape[1] == len(problem.masks)
    assert result.cross_success.shape[0] == len(problem.attacks) + len(result.additive_attacks)
    assert result.support_overlap.shape == (len(problem.masks), len(problem.masks))


def test_seeded_toy_case_is_deterministic():
    p1 = _toy_problem(seed=123)
    p2 = _toy_problem(seed=123)
    r1 = mpm.analyze_simultaneous_attacks(p1, include_additive=True, additive_max_terms=2)
    r2 = mpm.analyze_simultaneous_attacks(p2, include_additive=True, additive_max_terms=2)

    assert r1.cross_success.equals(r2.cross_success)
    assert np.allclose(r1.support_overlap.to_numpy(), r2.support_overlap.to_numpy())
    assert r1.dominance_edges == r2.dominance_edges


def test_dominance_edge_case_equal_coverage_prefers_smaller_norm():
    masks = [
        mpm.MaskDefense("∇1", np.array([1, 0], dtype=int), threshold=0.5),
        mpm.MaskDefense("∇2", np.array([0, 1], dtype=int), threshold=0.5),
    ]
    attacks = [
        mpm.AttackCandidate("Δ_small", np.array([0.6, 0.6], dtype=float)),
        mpm.AttackCandidate("Δ_large", np.array([1.3, 1.2], dtype=float)),
    ]
    problem = mpm.SimultaneousProblem(masks=masks, attacks=attacks)
    result = mpm.analyze_simultaneous_attacks(problem, include_additive=False)

    assert result.cross_success.loc["Δ_small"].tolist() == [1, 1]
    assert result.cross_success.loc["Δ_large"].tolist() == [1, 1]
    assert ("Δ_small", "Δ_large") in result.dominance_edges
    assert ("Δ_large", "Δ_small") not in result.dominance_edges
