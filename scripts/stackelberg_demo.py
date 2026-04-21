#!/usr/bin/env python3
"""Run a toy Stackelberg masked-game demo and export artifacts."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mpmgame import (
    StackelbergAttack,
    StackelbergInstance,
    StackelbergMask,
    attacker_best_response_to_mixed,
    export_stackelberg_outcome,
    solve_stackelberg_mixed,
    solve_stackelberg_pure,
)


def build_toy_instance() -> StackelbergInstance:
    masks = [
        StackelbergMask("∇0", cost=0.0),
        StackelbergMask("∇1", cost=1.0),
        StackelbergMask("∇2", cost=1.0),
        StackelbergMask("∇3", cost=2.0),
        StackelbergMask("∇4", cost=3.0),
    ]
    attacks = [StackelbergAttack(f"Δ{i+1}") for i in range(5)]

    # Attacker utility matrix U[defense, attack]. Defender minimizes these values.
    payoff = np.array(
        [
            [1.00, 0.95, 0.90, 1.10, 0.85],
            [0.60, 0.80, 0.65, 0.70, 0.75],
            [0.85, 0.55, 0.60, 0.72, 0.68],
            [0.45, 0.58, 0.52, 0.66, 0.61],
            [0.42, 0.40, 0.45, 0.59, 0.57],
        ],
        dtype=float,
    )
    return StackelbergInstance(payoff=payoff, masks=masks, attacks=attacks, budget=2.0)


def _plot_response_map(instance: StackelbergInstance, outpath: Path) -> None:
    best_attack_idx = np.argmax(instance.payoff, axis=1)
    best_attack_val = np.max(instance.payoff, axis=1)

    labels = [m.label for m in instance.masks]
    xticks = np.arange(len(labels))

    plt.figure(figsize=(7.2, 4.2))
    plt.scatter(xticks, best_attack_idx, c=best_attack_val, cmap="viridis", s=110)
    for i, (atk, val) in enumerate(zip(best_attack_idx, best_attack_val)):
        plt.text(i, atk + 0.08, f"{instance.attacks[atk].label}\n{val:.2f}", ha="center", va="bottom", fontsize=8)
    plt.xticks(xticks, labels)
    plt.yticks(np.arange(len(instance.attacks)), [a.label for a in instance.attacks])
    plt.xlabel("Defender committed pure mask")
    plt.ylabel("Attacker best response")
    plt.title("Pure-commitment response map")
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def _plot_value_comparison(mixed_value: float, pure_value: float, outpath: Path) -> None:
    plt.figure(figsize=(6.0, 4.0))
    names = ["Mixed commitment", "Pure commitment"]
    vals = [mixed_value, pure_value]
    colors = ["#1f77b4", "#ff7f0e"]
    plt.bar(names, vals, color=colors)
    plt.ylabel("Worst-case attacker value")
    plt.title("Stackelberg value comparison (lower is better)")
    for i, v in enumerate(vals):
        plt.text(i, v + 0.01, f"{v:.3f}", ha="center")
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def main() -> None:
    results_dir = Path("results/stackelberg")
    reports_dir = Path("reports/stackelberg")
    results_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    instance = build_toy_instance()
    mixed = solve_stackelberg_mixed(instance)
    pure = solve_stackelberg_pure(instance)

    export_stackelberg_outcome(instance, mixed, results_dir, stem="mixed")
    export_stackelberg_outcome(instance, pure, results_dir, stem="pure")

    best_idx, _, expected = attacker_best_response_to_mixed(instance.payoff, mixed.defender_strategy)
    expected_df = pd.DataFrame(
        {
            "attack_index": np.arange(len(instance.attacks)),
            "attack_label": [a.label for a in instance.attacks],
            "expected_attacker_utility": expected,
            "is_best_response": [i in set(best_idx) for i in range(len(instance.attacks))],
        }
    )
    expected_df.to_csv(results_dir / "mixed_expected_attack_utilities.csv", index=False)

    pd.DataFrame(instance.payoff, columns=[a.label for a in instance.attacks]).assign(
        defense_label=[m.label for m in instance.masks],
        defense_cost=[m.cost for m in instance.masks],
    ).to_csv(results_dir / "toy_payoff_matrix.csv", index=False)

    _plot_response_map(instance, reports_dir / "response_map.png")
    _plot_value_comparison(mixed.attacker_value, pure.attacker_value, reports_dir / "value_comparison.png")

    print(f"Saved Stackelberg demo outputs to {results_dir} and {reports_dir}")


if __name__ == "__main__":
    main()
