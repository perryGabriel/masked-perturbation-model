#!/usr/bin/env python3
"""Run a toy simultaneous-destabilization analysis and export artifacts."""

from __future__ import annotations

from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mpmgame import (
    AttackCandidate,
    MaskDefense,
    SimultaneousProblem,
    analyze_simultaneous_attacks,
    write_analysis_outputs,
)


def build_toy_problem(seed: int = 7) -> SimultaneousProblem:
    rng = np.random.default_rng(seed)
    dim = 6

    masks = [
        MaskDefense("∇1", np.array([1, 1, 0, 1, 0, 0], dtype=int), threshold=1.9),
        MaskDefense("∇2", np.array([1, 0, 1, 0, 1, 0], dtype=int), threshold=1.7),
        MaskDefense("∇3", np.array([0, 1, 1, 0, 0, 1], dtype=int), threshold=1.6),
        MaskDefense("∇4", np.array([1, 0, 0, 1, 0, 1], dtype=int), threshold=1.5),
    ]
    base = rng.uniform(0.2, 1.2, size=(4, dim))
    attacks = [AttackCandidate(f"Δ{i+1}", base[i]) for i in range(base.shape[0])]
    sensitivity = np.array([1.0, 1.2, 0.8, 1.1, 1.05, 0.95], dtype=float)
    return SimultaneousProblem(masks=masks, attacks=attacks, sensitivity=sensitivity)


def plot_matrix(mat: np.ndarray, xlabels: list[str], ylabels: list[str], title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    im = ax.imshow(mat, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(xlabels)), xlabels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(ylabels)), ylabels)
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_coverage_frequency(coverage: dict[str, list[str]], out_path: Path) -> None:
    labels = sorted(coverage)
    counts = [len(coverage[k]) for k in labels]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(labels, counts, color="#3C78D8")
    ax.set_ylabel("Number of masks destabilized")
    ax.set_xlabel("Attack candidate")
    ax.set_title("Destabilization frequency by attack")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> None:
    output_dir = Path("results/simultaneous_destabilization")
    output_dir.mkdir(parents=True, exist_ok=True)

    problem = build_toy_problem(seed=17)
    analysis = analyze_simultaneous_attacks(problem, include_additive=True, additive_max_terms=2)

    file_map = write_analysis_outputs(analysis, output_dir)

    plot_matrix(
        analysis.cross_success.to_numpy(dtype=float),
        list(analysis.cross_success.columns),
        list(analysis.cross_success.index),
        "Cross-success matrix (attack vs mask)",
        output_dir / "cross_success_heatmap.png",
    )
    plot_matrix(
        analysis.support_overlap.to_numpy(dtype=float),
        list(analysis.support_overlap.columns),
        list(analysis.support_overlap.index),
        "Support-overlap matrix (Jaccard)",
        output_dir / "support_overlap_heatmap.png",
    )
    plot_coverage_frequency(analysis.coverage, output_dir / "coverage_frequency.png")

    summary = {
        "num_masks": len(problem.masks),
        "num_base_attacks": len(problem.attacks),
        "num_additive_attacks": len(analysis.additive_attacks),
        "num_total_attacks": int(analysis.cross_success.shape[0]),
        "dominance_edges": analysis.dominance_edges,
        "output_files": {k: str(v) for k, v in file_map.items()},
        "figure_files": [
            str(output_dir / "cross_success_heatmap.png"),
            str(output_dir / "support_overlap_heatmap.png"),
            str(output_dir / "coverage_frequency.png"),
        ],
    }
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote simultaneous destabilization artifacts to: {output_dir}")


if __name__ == "__main__":
    main()
