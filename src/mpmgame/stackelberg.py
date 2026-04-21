"""Stackelberg defender-commitment utilities for finite masked games.

This module supports finite defender mask actions and finite attacker actions with
an attacker-utility payoff matrix ``U`` of shape ``(n_defenses, n_attacks)``.
The defender commits to a pure or mixed strategy over budget-feasible masks, and
the attacker best-responds by maximizing expected attacker utility.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.optimize import linprog


@dataclass(frozen=True)
class StackelbergMask:
    """Defender mask metadata for finite-action Stackelberg analysis."""

    label: str
    cost: float


@dataclass(frozen=True)
class StackelbergAttack:
    """Attacker action metadata for finite-action Stackelberg analysis."""

    label: str


@dataclass(frozen=True)
class StackelbergInstance:
    """Finite defender-commitment instance with optional budget constraint."""

    payoff: np.ndarray
    masks: list[StackelbergMask]
    attacks: list[StackelbergAttack]
    budget: float | None = None


@dataclass
class StackelbergOutcome:
    """Computed Stackelberg commitment outcome."""

    mode: str
    feasible_indices: list[int]
    defender_strategy: np.ndarray
    attacker_best_response_indices: list[int]
    attacker_value: float


def _validate_instance(instance: StackelbergInstance) -> None:
    if instance.payoff.ndim != 2:
        raise ValueError("payoff must be a 2D array")
    n_def, n_att = instance.payoff.shape
    if n_def == 0 or n_att == 0:
        raise ValueError("payoff must be non-empty")
    if len(instance.masks) != n_def:
        raise ValueError("number of masks must equal payoff rows")
    if len(instance.attacks) != n_att:
        raise ValueError("number of attacks must equal payoff columns")


def feasible_mask_indices(costs: np.ndarray, budget: float | None) -> list[int]:
    """Return defender-action indices satisfying the budget constraint."""

    if budget is None:
        return list(range(int(costs.size)))
    feasible = [int(i) for i, c in enumerate(costs.tolist()) if c <= float(budget) + 1e-12]
    if not feasible:
        raise ValueError("budget excludes all defender masks")
    return feasible


def attacker_best_response_to_mixed(
    payoff: np.ndarray,
    defender_mixed: np.ndarray,
    atol: float = 1e-10,
) -> tuple[list[int], float, np.ndarray]:
    """Compute attacker best responses to a defender mixed strategy."""

    expected = defender_mixed @ payoff
    value = float(np.max(expected))
    best = np.flatnonzero(expected >= value - atol).astype(int).tolist()
    return best, value, expected


def attacker_best_response_to_pure(payoff: np.ndarray, defender_index: int, atol: float = 1e-10) -> tuple[list[int], float]:
    """Compute attacker best responses to a defender pure commitment."""

    row = payoff[int(defender_index)]
    value = float(np.max(row))
    best = np.flatnonzero(row >= value - atol).astype(int).tolist()
    return best, value


def solve_stackelberg_mixed(instance: StackelbergInstance) -> StackelbergOutcome:
    """Solve defender mixed commitment under adversarial attacker tie-breaking."""

    _validate_instance(instance)
    payoff = instance.payoff.astype(float)
    costs = np.array([m.cost for m in instance.masks], dtype=float)
    feasible = feasible_mask_indices(costs, instance.budget)

    reduced = payoff[feasible, :]
    n_def, n_att = reduced.shape

    c = np.zeros(n_def + 1, dtype=float)
    c[-1] = 1.0

    a_ub = np.hstack([reduced.T, -np.ones((n_att, 1), dtype=float)])
    b_ub = np.zeros(n_att, dtype=float)

    a_eq = np.zeros((1, n_def + 1), dtype=float)
    a_eq[0, :n_def] = 1.0
    b_eq = np.array([1.0], dtype=float)

    bounds = [(0.0, None) for _ in range(n_def)] + [(None, None)]
    lp = linprog(c, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if not lp.success:
        raise RuntimeError(f"Stackelberg mixed LP failed: {lp.message}")

    strat_reduced = lp.x[:n_def]
    full_strat = np.zeros(payoff.shape[0], dtype=float)
    for k, i in enumerate(feasible):
        full_strat[i] = strat_reduced[k]

    best, value, _ = attacker_best_response_to_mixed(payoff, full_strat)
    return StackelbergOutcome(
        mode="mixed",
        feasible_indices=feasible,
        defender_strategy=full_strat,
        attacker_best_response_indices=best,
        attacker_value=value,
    )


def solve_stackelberg_pure(instance: StackelbergInstance) -> StackelbergOutcome:
    """Solve defender pure commitment under adversarial attacker tie-breaking."""

    _validate_instance(instance)
    payoff = instance.payoff.astype(float)
    costs = np.array([m.cost for m in instance.masks], dtype=float)
    feasible = feasible_mask_indices(costs, instance.budget)

    row_worst = payoff[feasible, :].max(axis=1)
    best_local = int(np.argmin(row_worst))
    defender_idx = feasible[best_local]

    defender_strategy = np.zeros(payoff.shape[0], dtype=float)
    defender_strategy[defender_idx] = 1.0
    best, value = attacker_best_response_to_pure(payoff, defender_idx)

    return StackelbergOutcome(
        mode="pure",
        feasible_indices=feasible,
        defender_strategy=defender_strategy,
        attacker_best_response_indices=best,
        attacker_value=value,
    )


def strategy_table(outcome: StackelbergOutcome, masks: list[StackelbergMask]) -> pd.DataFrame:
    """Convert defender commitment into a tabular policy representation."""

    rows = []
    feasible_set = set(outcome.feasible_indices)
    for i, m in enumerate(masks):
        rows.append(
            {
                "mask_index": i,
                "mask_label": m.label,
                "mask_cost": float(m.cost),
                "feasible": i in feasible_set,
                "probability": float(outcome.defender_strategy[i]),
            }
        )
    return pd.DataFrame(rows)


def export_stackelberg_outcome(
    instance: StackelbergInstance,
    outcome: StackelbergOutcome,
    output_dir: str | Path,
    stem: str,
) -> dict[str, Path]:
    """Write outcome summary and defender strategy table to disk."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    table_path = out / f"{stem}_policy_table.csv"
    summary_path = out / f"{stem}_summary.json"

    table = strategy_table(outcome, instance.masks)
    table.to_csv(table_path, index=False)

    payload = {
        "mode": outcome.mode,
        "budget": instance.budget,
        "feasible_indices": outcome.feasible_indices,
        "attacker_value": outcome.attacker_value,
        "attacker_best_response_indices": outcome.attacker_best_response_indices,
        "attacker_best_response_labels": [instance.attacks[j].label for j in outcome.attacker_best_response_indices],
        "defender_strategy": outcome.defender_strategy.tolist(),
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return {"policy_table": table_path, "summary": summary_path}
