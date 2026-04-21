import numpy as np

import mpmgame as mpm


def _success_sets() -> dict[str, set[str]]:
    return {
        "∇1": {"Δ1", "Δ2"},
        "∇2": {"Δ2", "Δ3"},
        "∇3": {"Δ1", "Δ2", "Δ3"},
        "∇4": {"Δ3"},
    }


def test_selection_size_respects_cap():
    w = mpm.MixedDefenseWeights(alpha_union=1.0, alpha_intersection=0.0, alpha_cardinality=0.3)
    sel = mpm.select_defense_subset_greedy(_success_sets(), z=2, weights=w)
    assert 1 <= len(sel.selected_labels) <= 2


def test_objective_monotonicity_under_zero_penalty_for_superset():
    sets = _success_sets()
    w = mpm.MixedDefenseWeights(alpha_union=1.0, alpha_intersection=0.5, alpha_cardinality=0.0)

    small = ["∇1"]
    large = ["∇1", "∇3"]
    assert set(sets["∇1"]).issubset(sets["∇3"])

    score_small = mpm.mixed_defense_objective(small, sets, w)
    score_large = mpm.mixed_defense_objective(large, sets, w)
    assert score_large >= score_small


def test_random_selection_seeded_behavior_is_deterministic():
    sets = _success_sets()
    w = mpm.MixedDefenseWeights(alpha_union=1.0, alpha_intersection=-0.25, alpha_cardinality=0.2)

    s1 = mpm.select_defense_subset_random(sets, z=3, weights=w, seed=42, num_trials=50)
    s2 = mpm.select_defense_subset_random(sets, z=3, weights=w, seed=42, num_trials=50)

    assert s1.selected_labels == s2.selected_labels
    assert np.isclose(s1.score, s2.score)


def test_reduced_game_evaluation_matches_selected_defenses():
    data = mpm.paper_example_data()
    sets = mpm.defense_success_sets(data.M, data.attacks, data.defenses)
    w = mpm.MixedDefenseWeights(alpha_union=1.0, alpha_intersection=0.0, alpha_cardinality=0.1)
    sel = mpm.select_defense_subset_greedy(sets, z=2, weights=w)

    eval_result = mpm.evaluate_mixed_defense_subset(data.M, data.attacks, data.defenses, sel)
    assert eval_result.reduced_payoff.shape[1] == len(sel.selected_labels)
