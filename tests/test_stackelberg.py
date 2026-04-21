import numpy as np

import mpmgame as mpm


def _toy_instance() -> mpm.StackelbergInstance:
    masks = [
        mpm.StackelbergMask("∇0", cost=0.0),
        mpm.StackelbergMask("∇1", cost=1.0),
        mpm.StackelbergMask("∇2", cost=2.0),
    ]
    attacks = [mpm.StackelbergAttack("Δ1"), mpm.StackelbergAttack("Δ2"), mpm.StackelbergAttack("Δ3")]
    payoff = np.array(
        [
            [0.9, 0.8, 1.0],
            [0.5, 0.7, 0.8],
            [0.6, 0.4, 0.6],
        ],
        dtype=float,
    )
    return mpm.StackelbergInstance(payoff=payoff, masks=masks, attacks=attacks, budget=1.0)


def test_best_response_to_mixed_matches_direct_max():
    instance = _toy_instance()
    strategy = np.array([0.2, 0.8, 0.0], dtype=float)

    best, value, expected = mpm.attacker_best_response_to_mixed(instance.payoff, strategy)

    assert np.isclose(value, expected.max())
    assert all(np.isclose(expected[idx], value) for idx in best)


def test_mixed_solution_objective_matches_best_response_value():
    instance = _toy_instance()
    outcome = mpm.solve_stackelberg_mixed(instance)

    feasible_set = set(outcome.feasible_indices)
    assert feasible_set == {0, 1}
    assert np.isclose(outcome.defender_strategy.sum(), 1.0)
    assert np.isclose(outcome.defender_strategy[2], 0.0)

    _, direct_value, _ = mpm.attacker_best_response_to_mixed(instance.payoff, outcome.defender_strategy)
    assert np.isclose(outcome.attacker_value, direct_value)


def test_mixed_commitment_is_no_worse_than_pure_commitment():
    instance = _toy_instance()
    mixed = mpm.solve_stackelberg_mixed(instance)
    pure = mpm.solve_stackelberg_pure(instance)

    assert mixed.attacker_value <= pure.attacker_value + 1e-10
