"""LFT/model-map construction utilities for MPM workflows.

The masked-perturbation game is posed in terms of a transfer map

    M(s) = C_r (sI - A)^{-1} B_w + D_rw,

where ``w`` is the write/attack channel and ``r`` is the read/monitor channel.
This module turns the state-space channel realization into the SISO-transfer-
function matrix representation used by :mod:`mpmgame.core`.
"""
from __future__ import annotations

import numpy as np
import control
from scipy.signal import ss2tf


TF = control.TransferFunction


def _as_2d_float(name: str, value: np.ndarray | list | tuple) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array")
    return arr


def _validate_model_map_shapes(
    A: np.ndarray | list | tuple,
    B_w: np.ndarray | list | tuple,
    C_r: np.ndarray | list | tuple,
    D_rw: np.ndarray | list | tuple,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    A_arr = _as_2d_float("A", A)
    B_arr = _as_2d_float("B_w", B_w)
    C_arr = _as_2d_float("C_r", C_r)
    D_arr = _as_2d_float("D_rw", D_rw)

    if A_arr.shape[0] != A_arr.shape[1]:
        raise ValueError("A must be square")
    n_x = A_arr.shape[0]
    if B_arr.shape[0] != n_x:
        raise ValueError("B_w row dimension must match A")
    if C_arr.shape[1] != n_x:
        raise ValueError("C_r column dimension must match A")
    if D_arr.shape != (C_arr.shape[0], B_arr.shape[1]):
        raise ValueError("D_rw must have shape (number of reads, number of writes)")

    return A_arr, B_arr, C_arr, D_arr


def build_model_map(
    A: np.ndarray | list | tuple,
    B_w: np.ndarray | list | tuple,
    C_r: np.ndarray | list | tuple,
    D_rw: np.ndarray | list | tuple,
) -> np.ndarray:
    """Construct the MPM model map ``M : w -> r`` as a transfer-function matrix.

    Parameters
    ----------
    A, B_w, C_r, D_rw:
        State-space matrices defining

        ``x_dot = A x + B_w w``
        ``r = C_r x + D_rw w``.

    Returns
    -------
    numpy.ndarray
        A 2D object array with shape ``(n_reads, n_writes)``. Each entry is a
        SISO :class:`control.TransferFunction`, so the result can be passed
        directly into the existing ``mpmgame`` attack/defense/game routines.

    Notes
    -----
    The construction is intentionally carried out one write channel at a time
    through :func:`scipy.signal.ss2tf`. This avoids relying on optional MIMO
    transfer-function conversion support while still producing the exact SISO
    entries of the LFT-derived map ``M(s)``.
    """
    A_arr, B_arr, C_arr, D_arr = _validate_model_map_shapes(A, B_w, C_r, D_rw)

    n_reads = C_arr.shape[0]
    n_writes = B_arr.shape[1]
    M = np.empty((n_reads, n_writes), dtype=object)

    for write_idx in range(n_writes):
        nums, den = ss2tf(A_arr, B_arr, C_arr, D_arr, input=write_idx)
        den = np.asarray(den, dtype=float)
        for read_idx in range(n_reads):
            num = np.asarray(nums[read_idx], dtype=float)
            M[read_idx, write_idx] = control.tf(num, den)

    return M


__all__ = ["build_model_map"]
