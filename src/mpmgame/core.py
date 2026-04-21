"""Core public API for masked-perturbation model games.

This module exposes functions centered on the paper-level objects:
- feedback map M
- attack operator Δ
- binary defense mask ∇
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Callable, Iterable, Sequence

import numpy as np

from .tf_tools import (
    tf_matrix,
    hadamard_tf,
    tf_add,
    tf_matmul,
    scale_rows_cols,
    is_unstable,
)
from .game import expected_utility as _expected_utility
from .game import solve_zero_sum_game as _solve_zero_sum_game

UtilityFn = Callable[[np.ndarray, np.ndarray], float]


@dataclass(frozen=True)
class AttackAction:
    label: str
    delta: np.ndarray


@dataclass(frozen=True)
class DefenseAction:
    label: str
    mask: np.ndarray


@dataclass
class SuccessSets:
    attack_success: dict[str, set[str]]
    defense_success: dict[str, set[str]]


@dataclass
class ReducedGame:
    payoff: np.ndarray
    attacks: list[AttackAction]
    defenses: list[DefenseAction]


@dataclass
class GameResult:
    full_payoff: np.ndarray
    reduced: ReducedGame
    attacker_mix: np.ndarray
    defender_mix: np.ndarray
    value: float
    dominated_attack_labels: set[str]
    dominated_defense_labels: set[str]
    success_sets: SuccessSets


def rank1_mask(eta_w: Sequence[int], eta_r: Sequence[int]) -> np.ndarray:
    eta_w_arr = np.asarray(eta_w, dtype=int).reshape(-1)
    eta_r_arr = np.asarray(eta_r, dtype=int).reshape(-1)
    return np.outer(eta_w_arr, eta_r_arr).astype(int)


def _recover_rank1_from_mask(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = np.asarray(mask, dtype=int)
    eta_w = (m.sum(axis=1) > 0).astype(int)
    eta_r = (m.sum(axis=0) > 0).astype(int)
    if not np.array_equal(rank1_mask(eta_w, eta_r), m):
        raise ValueError("mask is not rank-1 binary in η_w η_r^T form")
    return eta_w, eta_r


def defense_cost(
    mask: np.ndarray | None = None,
    eta_w: Sequence[int] | None = None,
    eta_r: Sequence[int] | None = None,
    c_w: Sequence[float] | None = None,
    c_r: Sequence[float] | None = None,
) -> float:
    if c_w is None or c_r is None:
        raise ValueError("c_w and c_r are required")
    cw = np.asarray(c_w, dtype=float)
    cr = np.asarray(c_r, dtype=float)
    if mask is not None:
        eta_w_arr, eta_r_arr = _recover_rank1_from_mask(np.asarray(mask, dtype=int))
    elif eta_w is not None and eta_r is not None:
        eta_w_arr = np.asarray(eta_w, dtype=int)
        eta_r_arr = np.asarray(eta_r, dtype=int)
    else:
        raise ValueError("provide either mask or (eta_w, eta_r)")
    return float((1 - eta_r_arr) @ cr + (1 - eta_w_arr) @ cw)


def enumerate_rank1_masks(num_writes: int, num_reads: int) -> list[np.ndarray]:
    masks = []
    for ew in product([0, 1], repeat=num_writes):
        for er in product([0, 1], repeat=num_reads):
            masks.append(rank1_mask(ew, er))
    return masks


def admissible_defenses(
    num_writes: int,
    num_reads: int,
    c_w: Sequence[float],
    c_r: Sequence[float],
    budget: float,
    include_empty: bool = False,
) -> list[np.ndarray]:
    out = []
    for m in enumerate_rank1_masks(num_writes, num_reads):
        if not include_empty and np.all(m == 0):
            continue
        if defense_cost(mask=m, c_w=c_w, c_r=c_r) <= budget:
            out.append(m)
    return out


def union_masks(mask1: np.ndarray, mask2: np.ndarray) -> np.ndarray:
    return (np.asarray(mask1, dtype=int) * np.asarray(mask2, dtype=int)).astype(int)


def complement_mask(mask: np.ndarray) -> np.ndarray:
    return (1 - np.asarray(mask, dtype=int)).astype(int)


def intersect_masks(mask1: np.ndarray, mask2: np.ndarray) -> np.ndarray:
    return complement_mask(union_masks(complement_mask(mask1), complement_mask(mask2)))


def is_subset_mask(mask_a: np.ndarray, mask_b: np.ndarray) -> bool:
    """True if defended set of a is subset of defended set of b."""
    da = complement_mask(mask_a)
    db = complement_mask(mask_b)
    return bool(np.all(da <= db))


def mask_footprint(mask: np.ndarray) -> dict[str, np.ndarray | int]:
    m = np.asarray(mask, dtype=int)
    defended = (m == 0)
    return {
        "defended_indices": np.argwhere(defended),
        "uncovered_indices": np.argwhere(m == 1),
        "num_defended": int(defended.sum()),
        "num_uncovered": int((m == 1).sum()),
    }


def hadamard_mask_delta(delta: np.ndarray, mask: np.ndarray) -> np.ndarray:
    dm = tf_matrix(delta)
    mm = np.asarray(mask, dtype=int)
    if dm.shape != mm.shape:
        raise ValueError("delta and mask must have same shape")
    out = np.empty(dm.shape, dtype=object)
    for i in range(dm.shape[0]):
        for j in range(dm.shape[1]):
            out[i, j] = int(mm[i, j]) * dm[i, j]
    return out


def masked_model_map(M: np.ndarray, mask_or_eta: np.ndarray | tuple[Sequence[int], Sequence[int]]) -> np.ndarray:
    Mm = tf_matrix(M)
    if isinstance(mask_or_eta, tuple):
        eta_w, eta_r = mask_or_eta
    else:
        eta_w, eta_r = _recover_rank1_from_mask(np.asarray(mask_or_eta, dtype=int))
    return scale_rows_cols(Mm, eta_r=eta_r, eta_w=eta_w)


def attack_sum(*deltas: np.ndarray) -> np.ndarray:
    return tf_add(*[tf_matrix(d) for d in deltas])


def closed_loop_interconnection(M: np.ndarray, delta: np.ndarray, mask: np.ndarray):
    Mm = tf_matrix(M)
    dmask = hadamard_mask_delta(delta, mask)
    return Mm, dmask


def is_destabilizing(
    M: np.ndarray,
    delta: np.ndarray,
    mask: np.ndarray,
    pole_tol: float = 1e-8,
    near_tol: float = 1e-6,
    method: str = "mask_delta",
) -> bool:
    Mm = tf_matrix(M)
    dmask = hadamard_mask_delta(delta, mask)
    if method == "mask_delta":
        return is_unstable(Mm, dmask, pole_tol=pole_tol, near_tol=near_tol)
    if method == "mask_model_map":
        Mmask = masked_model_map(Mm, mask)
        return is_unstable(Mmask, tf_matrix(delta), pole_tol=pole_tol, near_tol=near_tol)
    raise ValueError("method must be 'mask_delta' or 'mask_model_map'")


def utility(
    M: np.ndarray,
    delta: np.ndarray,
    mask: np.ndarray,
    ua: float = 1.0,
    ud: float = 0.0,
    ua_fn: UtilityFn | None = None,
    ud_fn: UtilityFn | None = None,
) -> float:
    if is_destabilizing(M, delta, mask):
        return float(ua_fn(delta, mask) if ua_fn is not None else ua)
    return float(ud_fn(delta, mask) if ud_fn is not None else ud)


def attack_success_set(M: np.ndarray, delta: np.ndarray, defenses: Sequence[np.ndarray]) -> list[int]:
    return [j for j, d in enumerate(defenses) if is_destabilizing(M, delta, d)]


def defense_success_set(M: np.ndarray, defense: np.ndarray, attacks: Sequence[np.ndarray]) -> list[int]:
    return [i for i, a in enumerate(attacks) if not is_destabilizing(M, a, defense)]


def compute_success_sets(
    M: np.ndarray,
    attacks: Sequence[AttackAction],
    defenses: Sequence[DefenseAction],
) -> SuccessSets:
    att = {
        a.label: {defenses[j].label for j in attack_success_set(M, a.delta, [d.mask for d in defenses])}
        for a in attacks
    }
    deff = {
        d.label: {attacks[i].label for i in defense_success_set(M, d.mask, [a.delta for a in attacks])}
        for d in defenses
    }
    return SuccessSets(attack_success=att, defense_success=deff)


def dominated_attacks(success_sets: dict[str, set[str]]) -> set[str]:
    labels = list(success_sets)
    dom: set[str] = set()
    for a in labels:
        for b in labels:
            if a == b:
                continue
            if success_sets[a] < success_sets[b]:
                dom.add(a)
    return dom


def dominated_defenses(success_sets: dict[str, set[str]]) -> set[str]:
    labels = list(success_sets)
    dom: set[str] = set()
    for a in labels:
        for b in labels:
            if a == b:
                continue
            if success_sets[a] < success_sets[b]:
                dom.add(a)
    return dom


def payoff_matrix(
    M: np.ndarray,
    attacks: Sequence[np.ndarray] | Sequence[AttackAction],
    defenses: Sequence[np.ndarray] | Sequence[DefenseAction],
    ua: float = 1.0,
    ud: float = 0.0,
    ua_fn: UtilityFn | None = None,
    ud_fn: UtilityFn | None = None,
) -> np.ndarray:
    attack_objs = [a.delta if isinstance(a, AttackAction) else a for a in attacks]
    defense_objs = [d.mask if isinstance(d, DefenseAction) else d for d in defenses]
    A = np.zeros((len(attack_objs), len(defense_objs)), dtype=float)
    for i, atk in enumerate(attack_objs):
        for j, dfn in enumerate(defense_objs):
            A[i, j] = utility(M, atk, dfn, ua=ua, ud=ud, ua_fn=ua_fn, ud_fn=ud_fn)
    return A


def eliminate_dominated_strategies(M: np.ndarray, attacks: Sequence[AttackAction], defenses: Sequence[DefenseAction]) -> ReducedGame:
    cur_att = list(attacks)
    cur_def = list(defenses)

    changed = True
    while changed:
        changed = False
        ss = compute_success_sets(M, cur_att, cur_def)
        d_att = dominated_attacks(ss.attack_success)
        d_def = dominated_defenses(ss.defense_success)
        if d_att:
            cur_att = [a for a in cur_att if a.label not in d_att]
            changed = True
        if d_def:
            cur_def = [d for d in cur_def if d.label not in d_def]
            changed = True

    reduced_payoff = payoff_matrix(M, cur_att, cur_def)
    return ReducedGame(payoff=reduced_payoff, attacks=cur_att, defenses=cur_def)




def select_defense_actions(defenses: Sequence[DefenseAction], labels: Sequence[str]) -> list[DefenseAction]:
    label_set = set(labels)
    selected = [d for d in defenses if d.label in label_set]
    if not selected:
        raise ValueError("no defenses selected")
    return selected


def evaluate_reduced_defense_game(
    M: np.ndarray,
    attacks: Sequence[np.ndarray] | Sequence[AttackAction],
    defenses: Sequence[np.ndarray] | Sequence[DefenseAction],
    ua: float = 1.0,
    ud: float = 0.0,
    ua_fn: UtilityFn | None = None,
    ud_fn: UtilityFn | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    reduced_payoff = payoff_matrix(M, attacks, defenses, ua=ua, ud=ud, ua_fn=ua_fn, ud_fn=ud_fn)
    attacker_mix, defender_mix, value = solve_zero_sum_game(reduced_payoff)
    return reduced_payoff, attacker_mix, defender_mix, value


def solve_zero_sum_game(payoff: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    return _solve_zero_sum_game(payoff)


def expected_utility(attacker_mix: np.ndarray, defender_mix: np.ndarray, payoff: np.ndarray) -> float:
    return _expected_utility(attacker_mix, defender_mix, payoff)
