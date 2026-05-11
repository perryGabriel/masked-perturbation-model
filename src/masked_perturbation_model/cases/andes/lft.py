"""LFT-style containers for linearized ANDES example systems.

This module intentionally provides only representation utilities.  It does not
construct simultaneous attack/destabilization operators.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class AndesLFT:
    """State-space/LFT-style representation of a linearized ANDES case.

    The representation records compatible state-space matrices and identity
    read/write masks over the linearized states.  It is deliberately a data
    container for larger examples, not a synthesis routine.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    read_mask: np.ndarray
    write_mask: np.ndarray
    state_names: tuple[str, ...] = ()
    input_names: tuple[str, ...] = ()
    output_names: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def nstates(self) -> int:
        return int(self.A.shape[0])

    @property
    def ninputs(self) -> int:
        return int(self.B.shape[1])

    @property
    def noutputs(self) -> int:
        return int(self.C.shape[0])

    @property
    def mask_shape(self) -> tuple[int, int]:
        return (int(self.write_mask.shape[0]), int(self.read_mask.shape[0]))


def build_lft_from_state_matrix(
    A: np.ndarray,
    *,
    state_names: tuple[str, ...] = (),
    metadata: Mapping[str, Any] | None = None,
) -> AndesLFT:
    """Build a deterministic LFT-style container from a square state matrix."""

    A_arr = np.asarray(A, dtype=float)
    if A_arr.ndim != 2 or A_arr.shape[0] != A_arr.shape[1]:
        raise ValueError("A must be a square two-dimensional state matrix")
    n = int(A_arr.shape[0])
    names = state_names or tuple(f"x{i}" for i in range(n))
    if len(names) != n:
        raise ValueError("state_names length must match A dimensions")
    B = np.eye(n)
    C = np.eye(n)
    D = np.zeros((n, n))
    return AndesLFT(
        A=A_arr,
        B=B,
        C=C,
        D=D,
        read_mask=np.ones(n, dtype=int),
        write_mask=np.ones(n, dtype=int),
        state_names=tuple(names),
        input_names=tuple(f"u_{name}" for name in names),
        output_names=tuple(f"y_{name}" for name in names),
        metadata=dict(metadata or {}),
    )
