#!/usr/bin/env python3
"""Benchmark mixed-defense action-set selection over z and objective weights."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mpmgame import (
    MixedDefenseWeights,
    defense_success_sets,
    evaluate_mixed_defense_subset,
    paper_example_data,
    payoff_matrix,
    select_defense_subset_greedy,
    select_defense_subset_random,
    solve_zero_sum_game,
)


def main() -> None:
    results_dir = Path("results/mixed_defense")
    reports_dir = Path("reports/mixed_defense")
    results_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    data = paper_example_data()
    success_sets = defense_success_sets(data.M, data.attacks, data.defenses)

    full_pay = payoff_matrix(data.M, data.attacks, data.defenses)
    _, _, full_value = solve_zero_sum_game(full_pay)

    z_values = [1, 2, 3, 4]
    weight_grid = [
        MixedDefenseWeights(alpha_union=1.0, alpha_intersection=0.0, alpha_cardinality=0.0),
        MixedDefenseWeights(alpha_union=1.0, alpha_intersection=0.5, alpha_cardinality=0.2),
        MixedDefenseWeights(alpha_union=1.0, alpha_intersection=-0.25, alpha_cardinality=0.4),
    ]

    rows: list[dict[str, object]] = []
    for z in z_values:
        for idx, weights in enumerate(weight_grid):
            for method in ("greedy", "random"):
                if method == "greedy":
                    selection = select_defense_subset_greedy(success_sets, z=z, weights=weights)
                else:
                    selection = select_defense_subset_random(
                        success_sets,
                        z=z,
                        weights=weights,
                        seed=100 + z * 10 + idx,
                        num_trials=400,
                    )
                eval_result = evaluate_mixed_defense_subset(data.M, data.attacks, data.defenses, selection)
                union_coverage = len(set().union(*[success_sets[k] for k in selection.selected_labels]))
                rows.append(
                    {
                        "method": method,
                        "z": z,
                        "weights_id": idx,
                        "alpha_union": weights.alpha_union,
                        "alpha_intersection": weights.alpha_intersection,
                        "alpha_cardinality": weights.alpha_cardinality,
                        "selected_labels": ",".join(selection.selected_labels),
                        "selected_size": len(selection.selected_labels),
                        "objective_score": selection.score,
                        "pair_union": selection.components.pair_union,
                        "pair_intersection": selection.components.pair_intersection,
                        "cardinality_penalty": selection.components.cardinality_penalty,
                        "union_coverage": union_coverage,
                        "reduced_value": eval_result.value,
                        "full_value": full_value,
                        "value_gap": eval_result.value - full_value,
                    }
                )

    df = pd.DataFrame(rows).sort_values(["weights_id", "z", "method"])
    df.to_csv(results_dir / "sweep_results.csv", index=False)

    summary = {
        "full_game_value": full_value,
        "num_attacks": len(data.attacks),
        "num_defenses": len(data.defenses),
        "z_values": z_values,
        "weights": [asdict(w) for w in weight_grid],
        "best_by_method": {
            method: df[df["method"] == method].sort_values("reduced_value", ascending=False).head(1).to_dict("records")[0]
            for method in ("greedy", "random")
        },
    }
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    pivot_val = df.pivot_table(index="z", columns="method", values="reduced_value", aggfunc="mean")
    ax = pivot_val.plot(marker="o", figsize=(7, 4.5), title="Reduced-game value vs z")
    ax.axhline(full_value, linestyle="--", color="k", linewidth=1.0, label="full-game value")
    ax.set_ylabel("Game value")
    ax.legend()
    plt.tight_layout()
    plt.savefig(reports_dir / "value_vs_z.png", dpi=160)
    plt.close()

    pivot_cov = df.pivot_table(index="z", columns="method", values="union_coverage", aggfunc="mean")
    ax = pivot_cov.plot(marker="s", figsize=(7, 4.5), title="Union coverage vs z")
    ax.set_ylabel("Covered attacks")
    plt.tight_layout()
    plt.savefig(reports_dir / "coverage_vs_z.png", dpi=160)
    plt.close()

    ax = df.plot.scatter(
        x="objective_score",
        y="reduced_value",
        c="z",
        cmap="viridis",
        figsize=(7, 4.5),
        title="Objective score vs reduced-game value",
    )
    ax.axhline(full_value, linestyle="--", color="k", linewidth=1.0)
    plt.tight_layout()
    plt.savefig(reports_dir / "objective_vs_value_scatter.png", dpi=160)
    plt.close()

    print(f"Saved benchmark outputs to {results_dir} and plots to {reports_dir}")


if __name__ == "__main__":
    main()
