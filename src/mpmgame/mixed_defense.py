"""Mixed-defense action-set design via success-set geometry."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import numpy as np

from .core import AttackAction, DefenseAction, compute_success_sets, evaluate_reduced_defense_game


@dataclass(frozen=True)
class MixedDefenseWeights:
    alpha_union: float = 1.0
    alpha_intersection: float = 0.0
    alpha_cardinality: float = 0.0


@dataclass(frozen=True)
class MixedDefenseObjectiveComponents:
    pair_union: float
    pair_intersection: float
    cardinality_penalty: float


@dataclass
class MixedDefenseSelection:
    method: str
    selected_labels: list[str]
    score: float
    components: MixedDefenseObjectiveComponents


@dataclass
class MixedDefenseEvaluation:
    selection: MixedDefenseSelection
    reduced_payoff: np.ndarray
    attacker_mix: np.ndarray
    defender_mix: np.ndarray
    value: float


def defense_success_sets(
    M: np.ndarray,
    attacks: Sequence[AttackAction],
    defenses: Sequence[DefenseAction],
) -> dict[str, set[str]]:
    """Return S(∇): attack labels each defense successfully blocks."""
    return compute_success_sets(M, attacks, defenses).defense_success


def _pair_scores(selected_labels: Sequence[str], success_sets: dict[str, set[str]]) -> tuple[float, float]:
    labels = list(selected_labels)
    if len(labels) == 1:
        s = success_sets[labels[0]]
        return float(len(s)), float(len(s))

    pair_union = 0.0
    pair_intersection = 0.0
    for a, b in combinations(labels, 2):
        sa = success_sets[a]
        sb = success_sets[b]
        pair_union += len(sa | sb)
        pair_intersection += len(sa & sb)
    return pair_union, pair_intersection


def mixed_defense_objective_components(
    selected_labels: Sequence[str],
    success_sets: dict[str, set[str]],
    weights: MixedDefenseWeights,
) -> MixedDefenseObjectiveComponents:
    pair_union, pair_intersection = _pair_scores(selected_labels, success_sets)
    return MixedDefenseObjectiveComponents(
        pair_union=pair_union,
        pair_intersection=pair_intersection,
        cardinality_penalty=weights.alpha_cardinality * len(selected_labels),
    )


def mixed_defense_objective(
    selected_labels: Sequence[str],
    success_sets: dict[str, set[str]],
    weights: MixedDefenseWeights,
) -> float:
    if not selected_labels:
        return -np.inf
    c = mixed_defense_objective_components(selected_labels, success_sets, weights)
    return float(weights.alpha_union * c.pair_union + weights.alpha_intersection * c.pair_intersection - c.cardinality_penalty)


def select_defense_subset_greedy(
    success_sets: dict[str, set[str]],
    z: int,
    weights: MixedDefenseWeights,
) -> MixedDefenseSelection:
    if z <= 0:
        raise ValueError("z must be positive")

    available = set(success_sets)
    selected: list[str] = []

    while available and len(selected) < z:
        cur_score = mixed_defense_objective(selected, success_sets, weights) if selected else -np.inf
        best_label = None
        best_score = cur_score
        for label in sorted(available):
            cand = selected + [label]
            cand_score = mixed_defense_objective(cand, success_sets, weights)
            if cand_score > best_score:
                best_score = cand_score
                best_label = label
        if best_label is None:
            break
        selected.append(best_label)
        available.remove(best_label)

    if not selected:
        selected = [sorted(success_sets)[0]]

    components = mixed_defense_objective_components(selected, success_sets, weights)
    score = mixed_defense_objective(selected, success_sets, weights)
    return MixedDefenseSelection("greedy", selected, score, components)


def select_defense_subset_random(
    success_sets: dict[str, set[str]],
    z: int,
    weights: MixedDefenseWeights,
    seed: int = 0,
    num_trials: int = 100,
) -> MixedDefenseSelection:
    if z <= 0:
        raise ValueError("z must be positive")
    labels = sorted(success_sets)
    rng = np.random.default_rng(seed)

    best_labels = labels[:1]
    best_score = mixed_defense_objective(best_labels, success_sets, weights)

    max_k = min(z, len(labels))
    for _ in range(num_trials):
        k = int(rng.integers(1, max_k + 1))
        draw = sorted(rng.choice(labels, size=k, replace=False).tolist())
        score = mixed_defense_objective(draw, success_sets, weights)
        if score > best_score:
            best_score = score
            best_labels = draw

    components = mixed_defense_objective_components(best_labels, success_sets, weights)
    return MixedDefenseSelection("random", best_labels, best_score, components)


def evaluate_mixed_defense_subset(
    M: np.ndarray,
    attacks: Sequence[AttackAction],
    defenses: Sequence[DefenseAction],
    selection: MixedDefenseSelection,
) -> MixedDefenseEvaluation:
    reduced_defenses = [d for d in defenses if d.label in set(selection.selected_labels)]
    reduced_payoff, attacker_mix, defender_mix, value = evaluate_reduced_defense_game(
        M=M,
        attacks=attacks,
        defenses=reduced_defenses,
    )
    return MixedDefenseEvaluation(
        selection=selection,
        reduced_payoff=reduced_payoff,
        attacker_mix=attacker_mix,
        defender_mix=defender_mix,
        value=value,
    )
