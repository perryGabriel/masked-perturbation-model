"""Small-signal extraction and fallback systems."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class LinearModel:
    A: np.ndarray
    state_names: list[str]
    source: str


def verify_stability(A: np.ndarray, tol: float = 1e-8) -> tuple[bool, np.ndarray]:
    """Return strict Hurwitz stability status and eigenvalues of ``A``."""
    eigvals = np.linalg.eigvals(A)
    stable = np.max(np.real(eigvals)) < -tol
    return stable, eigvals


def extract_andes_state_matrix(system: Any) -> LinearModel:
    """Try to extract a reduced small-signal state matrix from ANDES.

    This remains intentionally defensive because exact attributes can vary with
    ANDES version and case format. The demo notebook surfaces diagnostics if no
    state matrix can be located.
    """
    if hasattr(system, "EIG"):
        system.EIG.run()

    candidate_objects = []
    if hasattr(system, "EIG"):
        candidate_objects.append(system.EIG)
    candidate_objects.append(system)

    matrix_names = ["As", "A", "state_matrix", "Amat", "A_matrix"]
    name_names = ["state_names", "x_names", "states", "x_name"]

    A = None
    for obj in candidate_objects:
        for name in matrix_names:
            if hasattr(obj, name):
                value = getattr(obj, name)
                try:
                    arr = np.asarray(value, dtype=float)
                except Exception:
                    continue
                if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] > 0:
                    A = arr
                    break
        if A is not None:
            break

    if A is None:
        raise AttributeError(
            "Could not locate an ANDES small-signal state matrix. "
            "Inspect system.EIG and patch extract_andes_state_matrix()."
        )

    state_names: list[str] | None = None
    for obj in candidate_objects:
        for name in name_names:
            if hasattr(obj, name):
                value = getattr(obj, name)
                try:
                    state_names = [str(v) for v in list(value)]
                except Exception:
                    state_names = None
                if state_names and len(state_names) == A.shape[0]:
                    break
        if state_names and len(state_names) == A.shape[0]:
            break

    if state_names is None or len(state_names) != A.shape[0]:
        state_names = [f"x{k}" for k in range(A.shape[0])]

    return LinearModel(A=A, state_names=state_names, source="ANDES")


def synthetic_two_area_like_model() -> LinearModel:
    """Create a stable four-machine swing-like fallback model.

    This is not a benchmark substitute. It keeps the MPM verification path
    executable while case-specific ANDES extraction details are refined.
    """
    n_gen = 4
    M = np.array([4.0, 4.5, 4.0, 4.2])
    D = np.array([1.2, 1.1, 1.3, 1.0])

    L = np.array(
        [
            [6.0, -3.0, -2.0, -1.0],
            [-3.0, 6.0, -1.0, -2.0],
            [-2.0, -1.0, 6.0, -3.0],
            [-1.0, -2.0, -3.0, 6.0],
        ]
    )

    K_delta = 0.4 * np.eye(n_gen)
    K_omega = 0.8 * np.eye(n_gen)

    Z = np.zeros((n_gen, n_gen))
    identity = np.eye(n_gen)
    Minv = np.diag(1.0 / M)

    A = np.block(
        [
            [Z, identity],
            [-Minv @ (L + K_delta), -Minv @ (np.diag(D) + K_omega)],
        ]
    )

    A[:n_gen, :n_gen] -= 0.02 * np.eye(n_gen)

    state_names = [f"delta_g{k+1}" for k in range(n_gen)] + [f"omega_g{k+1}" for k in range(n_gen)]
    return LinearModel(A=A, state_names=state_names, source="synthetic_fallback")
