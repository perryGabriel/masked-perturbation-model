"""Simultaneous destabilization utilities for mask sets and candidate attacks.

This module provides a lightweight exploratory workflow focused on finite attack
candidates and finite defense masks. A per-mask destabilization model can be
plugged in as a callable, and a default linear-threshold model is provided for
reproducible toy experiments.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import json

import numpy as np
import pandas as pd


ArrayLike1D = np.ndarray


@dataclass(frozen=True)
class MaskDefense:
    """Single defense mask over an abstract attack-support basis."""

    label: str
    mask: ArrayLike1D
    threshold: float = 1.0


@dataclass(frozen=True)
class AttackCandidate:
    """Single attack candidate with optional support metadata."""

    label: str
    vector: ArrayLike1D


@dataclass(frozen=True)
class SimultaneousProblem:
    """Finite simultaneous-destabilization problem definition."""

    masks: list[MaskDefense]
    attacks: list[AttackCandidate]
    sensitivity: ArrayLike1D | None = None


@dataclass
class SimultaneousAnalysisResult:
    """Computed matrices and relations for a simultaneous attack study."""

    cross_success: pd.DataFrame
    support_overlap: pd.DataFrame
    dominance_edges: list[tuple[str, str]]
    coverage: dict[str, list[str]]
    additive_attacks: list[AttackCandidate]


def _validate_lengths(problem: SimultaneousProblem) -> int:
    if not problem.masks:
        raise ValueError("problem.masks must be non-empty")
    if not problem.attacks:
        raise ValueError("problem.attacks must be non-empty")
    d = int(problem.masks[0].mask.size)
    for m in problem.masks:
        if m.mask.size != d:
            raise ValueError("all mask vectors must have identical dimension")
    for a in problem.attacks:
        if a.vector.size != d:
            raise ValueError("all attack vectors must match mask dimension")
    if problem.sensitivity is not None and problem.sensitivity.size != d:
        raise ValueError("sensitivity vector must match mask dimension")
    return d


def destabilizes(mask: MaskDefense, attack: AttackCandidate, sensitivity: ArrayLike1D | None = None) -> bool:
    """Evaluate default per-mask destabilization success for one attack."""

    s = np.ones_like(mask.mask, dtype=float) if sensitivity is None else sensitivity.astype(float)
    score = float(np.dot(mask.mask.astype(float) * s, np.abs(attack.vector.astype(float))))
    return score >= float(mask.threshold)


def evaluate_attack_against_masks(
    attack: AttackCandidate,
    masks: list[MaskDefense],
    sensitivity: ArrayLike1D | None = None,
) -> dict[str, bool]:
    """Return per-mask destabilization outcomes for one attack."""

    return {m.label: destabilizes(m, attack, sensitivity=sensitivity) for m in masks}


def cross_success_matrix(problem: SimultaneousProblem) -> pd.DataFrame:
    """Build attack-by-mask boolean success matrix."""

    _validate_lengths(problem)
    rows = []
    for attack in problem.attacks:
        outcome = evaluate_attack_against_masks(attack, problem.masks, sensitivity=problem.sensitivity)
        rows.append([int(outcome[m.label]) for m in problem.masks])
    return pd.DataFrame(rows, index=[a.label for a in problem.attacks], columns=[m.label for m in problem.masks])


def support_overlap_matrix(masks: list[MaskDefense], normalize: bool = True) -> pd.DataFrame:
    """Build pairwise mask-support overlap matrix (Jaccard by default)."""

    labels = [m.label for m in masks]
    arr = np.zeros((len(masks), len(masks)), dtype=float)
    for i, mi in enumerate(masks):
        si = mi.mask.astype(bool)
        for j, mj in enumerate(masks):
            sj = mj.mask.astype(bool)
            inter = np.logical_and(si, sj).sum()
            if not normalize:
                arr[i, j] = float(inter)
                continue
            union = np.logical_or(si, sj).sum()
            arr[i, j] = 1.0 if union == 0 else float(inter) / float(union)
    return pd.DataFrame(arr, index=labels, columns=labels)


def dominance_relations(
    cross_success: pd.DataFrame,
    attack_vectors: dict[str, ArrayLike1D] | None = None,
    tol: float = 1e-12,
) -> list[tuple[str, str]]:
    """Return dominance edges (dominator, dominated) over attack candidates.

    An attack A dominates B if success(A) is a superset of success(B). If both
    are equal, tie-breaking prefers lower or equal L1 norm when vectors exist.
    """

    labels = list(cross_success.index)
    edges: list[tuple[str, str]] = []
    for i, a in enumerate(labels):
        va = cross_success.loc[a].to_numpy(dtype=int)
        for j, b in enumerate(labels):
            if i == j:
                continue
            vb = cross_success.loc[b].to_numpy(dtype=int)
            strict_superset = np.all(va >= vb) and np.any(va > vb)
            equal_rows = np.array_equal(va, vb)
            if strict_superset:
                edges.append((a, b))
                continue
            if equal_rows and attack_vectors is not None:
                norm_a = np.linalg.norm(attack_vectors[a], ord=1)
                norm_b = np.linalg.norm(attack_vectors[b], ord=1)
                if norm_a <= norm_b + tol and (norm_a < norm_b - tol or a < b):
                    edges.append((a, b))
    return sorted(set(edges))


def additive_combination_search(
    attacks: list[AttackCandidate],
    max_terms: int = 2,
    max_combinations: int | None = None,
) -> list[AttackCandidate]:
    """Generate optional additive attack combinations from the given candidates."""

    if max_terms < 2:
        return []
    generated: list[AttackCandidate] = []
    for r in range(2, max_terms + 1):
        for combo in combinations(attacks, r):
            label = " + ".join(a.label for a in combo)
            vec = np.sum([a.vector for a in combo], axis=0)
            generated.append(AttackCandidate(label=label, vector=vec))
            if max_combinations is not None and len(generated) >= max_combinations:
                return generated
    return generated


def analyze_simultaneous_attacks(
    problem: SimultaneousProblem,
    include_additive: bool = False,
    additive_max_terms: int = 2,
) -> SimultaneousAnalysisResult:
    """Run end-to-end simultaneous attack analysis for one finite problem."""

    _validate_lengths(problem)
    additive = additive_combination_search(problem.attacks, max_terms=additive_max_terms) if include_additive else []
    all_attacks = [*problem.attacks, *additive]
    expanded = SimultaneousProblem(masks=problem.masks, attacks=all_attacks, sensitivity=problem.sensitivity)
    cross = cross_success_matrix(expanded)
    overlap = support_overlap_matrix(problem.masks)
    vec_map = {a.label: a.vector for a in all_attacks}
    edges = dominance_relations(cross, attack_vectors=vec_map)
    coverage = {
        label: [m for m, v in row.items() if int(v) == 1]
        for label, row in cross.to_dict(orient="index").items()
    }
    return SimultaneousAnalysisResult(
        cross_success=cross,
        support_overlap=overlap,
        dominance_edges=edges,
        coverage=coverage,
        additive_attacks=additive,
    )


def write_analysis_outputs(result: SimultaneousAnalysisResult, output_dir: str | Path) -> dict[str, Path]:
    """Write machine-readable analysis outputs (CSV/JSON)."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    cross_csv = out / "cross_success_matrix.csv"
    overlap_csv = out / "support_overlap_matrix.csv"
    dom_json = out / "dominance_relations.json"
    cov_json = out / "coverage_by_attack.json"
    additive_json = out / "additive_attacks.json"

    result.cross_success.to_csv(cross_csv)
    result.support_overlap.to_csv(overlap_csv)
    dom_payload = [{"dominant": a, "dominated": b} for a, b in result.dominance_edges]
    cov_payload = result.coverage
    additive_payload = [{"label": a.label, "vector": a.vector.tolist()} for a in result.additive_attacks]

    dom_json.write_text(json.dumps(dom_payload, indent=2), encoding="utf-8")
    cov_json.write_text(json.dumps(cov_payload, indent=2), encoding="utf-8")
    additive_json.write_text(json.dumps(additive_payload, indent=2), encoding="utf-8")

    return {
        "cross_success_csv": cross_csv,
        "support_overlap_csv": overlap_csv,
        "dominance_json": dom_json,
        "coverage_json": cov_json,
        "additive_json": additive_json,
    }
