"""Paper example construction and reproduction helpers."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import control

from .core import (
    AttackAction,
    DefenseAction,
    GameResult,
    SuccessSets,
    admissible_defenses,
    attack_sum,
    compute_success_sets,
    dominated_attacks,
    dominated_defenses,
    eliminate_dominated_strategies,
    is_destabilizing,
    payoff_matrix,
    solve_zero_sum_game,
)
from .tf_tools import tf_matrix


@dataclass
class PaperExampleData:
    M: np.ndarray
    attacks: list[AttackAction]
    defenses: list[DefenseAction]
    c_r: np.ndarray
    c_w: np.ndarray
    budget: float


def paper_example_data() -> PaperExampleData:
    s = control.tf([1, 0], [1])
    P = 3 / (3 * s + 1)
    K = 3 / (3 * s + 10)
    den = 1 + 8 * P * K

    M = tf_matrix([
        [P / den, (P * K) / den],
        [(-8 * P * K) / den, K / den],
    ])

    d1 = tf_matrix([
        [(2.8 * s - 1.62) / (s + 0.58), 0],
        [0, 0],
    ])
    d2 = tf_matrix([
        [0, 0],
        [0, (4.7 * s - 0.8) / (s + 0.001)],
    ])
    d3 = tf_matrix([
        [0, (-1.2 * s + 0.7) / (s + 0.6)],
        [0, 0],
    ])

    attacks = [
        AttackAction("Δ1", d1),
        AttackAction("Δ2", d2),
        AttackAction("Δ3", d3),
    ]

    n1 = np.array([[1, 0], [1, 0]], dtype=int)
    n2 = np.array([[1, 1], [0, 0]], dtype=int)
    n3 = np.array([[0, 0], [1, 1]], dtype=int)
    n4 = np.array([[1, 0], [0, 0]], dtype=int)
    defenses = [
        DefenseAction("∇1", n1),
        DefenseAction("∇2", n2),
        DefenseAction("∇3", n3),
        DefenseAction("∇4", n4),
    ]

    c_r = np.array([3, 1], dtype=float)
    c_w = np.array([2, 1], dtype=float)
    budget = 2.0

    return PaperExampleData(M=M, attacks=attacks, defenses=defenses, c_r=c_r, c_w=c_w, budget=budget)


def run_paper_example() -> GameResult:
    data = paper_example_data()
    success = compute_success_sets(data.M, data.attacks, data.defenses)
    full_payoff = payoff_matrix(data.M, data.attacks, data.defenses, ua=1.0, ud=0.0)

    dom_att = dominated_attacks(success.attack_success)
    dom_def = dominated_defenses(success.defense_success)

    reduced = eliminate_dominated_strategies(data.M, data.attacks, data.defenses)
    att_mix, def_mix, value = solve_zero_sum_game(reduced.payoff)

    return GameResult(
        full_payoff=full_payoff,
        reduced=reduced,
        attacker_mix=att_mix,
        defender_mix=def_mix,
        value=value,
        dominated_attack_labels=dom_att,
        dominated_defense_labels=dom_def,
        success_sets=success,
    )


def coordinated_attack_defeats_all() -> bool:
    data = paper_example_data()
    d4 = attack_sum(data.attacks[0].delta, data.attacks[1].delta)
    return all(is_destabilizing(data.M, d4, d.mask) for d in data.defenses)
