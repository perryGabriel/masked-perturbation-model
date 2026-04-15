"""Transfer-function helpers for masked-perturbation games.

The public conceptual layer remains in transfer maps M, attacks Δ, and masks ∇.
"""

from __future__ import annotations

from typing import Iterable, Sequence
import warnings

import numpy as np
import control
from control.exception import ControlMIMONotImplemented


TF = control.TransferFunction


def _as_tf(x: object) -> TF:
    if isinstance(x, control.TransferFunction):
        return x
    if isinstance(x, (int, float, complex, np.number)):
        return control.tf([x], [1])
    return control.tf(x)


def tf_matrix(entries: Sequence[Sequence[object]]) -> np.ndarray:
    """Create a 2D ndarray of SISO transfer functions."""
    rows = []
    for row in entries:
        rows.append([_as_tf(v) for v in row])
    arr = np.array(rows, dtype=object)
    if arr.ndim != 2:
        raise ValueError("entries must define a 2D matrix")
    return arr


def zeros_tf(shape: tuple[int, int]) -> np.ndarray:
    return np.array([[control.tf([0], [1]) for _ in range(shape[1])] for _ in range(shape[0])], dtype=object)


def eye_tf(n: int) -> np.ndarray:
    out = zeros_tf((n, n))
    one = control.tf([1], [1])
    for i in range(n):
        out[i, i] = one
    return out


def hadamard_tf(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    if A.shape != B.shape:
        raise ValueError("Hadamard operands must have the same shape")
    out = zeros_tf(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            out[i, j] = _as_tf(A[i, j]) * _as_tf(B[i, j])
    return out


def tf_add(*mats: np.ndarray) -> np.ndarray:
    if not mats:
        raise ValueError("at least one matrix required")
    out = tf_matrix(mats[0])
    for mat in mats[1:]:
        if out.shape != mat.shape:
            raise ValueError("all matrices must have same shape")
        for i in range(out.shape[0]):
            for j in range(out.shape[1]):
                out[i, j] = out[i, j] + _as_tf(mat[i, j])
    return out


def tf_matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    if A.shape[1] != B.shape[0]:
        raise ValueError("matrix dimensions do not align")
    out = zeros_tf((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            acc = control.tf([0], [1])
            for k in range(A.shape[1]):
                acc = acc + _as_tf(A[i, k]) * _as_tf(B[k, j])
            out[i, j] = acc
    return out


def scale_rows_cols(M: np.ndarray, eta_r: Sequence[int], eta_w: Sequence[int]) -> np.ndarray:
    """Construct M_∇ = diag(η_r) M diag(η_w)."""
    eta_r_arr = np.asarray(eta_r, dtype=int).reshape(-1)
    eta_w_arr = np.asarray(eta_w, dtype=int).reshape(-1)
    if M.shape != (eta_r_arr.size, eta_w_arr.size):
        raise ValueError("eta dimensions must match M shape")
    out = zeros_tf(M.shape)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            out[i, j] = int(eta_r_arr[i]) * _as_tf(M[i, j]) * int(eta_w_arr[j])
    return out


def to_mimo_tf(M: np.ndarray) -> TF:
    """Convert object-array transfer matrix to python-control MIMO TransferFunction."""
    p, q = M.shape
    num = [[None for _ in range(q)] for _ in range(p)]
    den = [[None for _ in range(q)] for _ in range(p)]
    for i in range(p):
        for j in range(q):
            gij = _as_tf(M[i, j])
            num[i][j] = np.array(gij.num[0][0], dtype=float)
            den[i][j] = np.array(gij.den[0][0], dtype=float)
    return control.TransferFunction(num, den)


def feedback_resolvent(M: np.ndarray, Delta: np.ndarray) -> control.StateSpace:
    """Return (I - MΔ)^(-1) as a dynamic system for stability checks."""
    if M.shape[1] != Delta.shape[0] or M.shape[0] != Delta.shape[1]:
        raise ValueError("Expected M shape p×q and Delta shape q×p")
    L = tf_matmul(M, Delta)
    I = eye_tf(L.shape[0])
    # `python-control` does not implement MIMO feedback for TransferFunction
    # objects in some versions, but the same interconnection works for
    # StateSpace systems. Convert before closing the loop.
    I_sys = control.ss(to_mimo_tf(I))
    L_sys = control.ss(to_mimo_tf(L))
    return control.feedback(I_sys, L_sys, sign=1)


def _tf_det(A: np.ndarray) -> TF:
    """Determinant of a square transfer-matrix as a SISO transfer function."""
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("determinant requires a square 2D matrix")
    n = A.shape[0]
    if n == 1:
        return _as_tf(A[0, 0])
    out = control.tf([0], [1])
    for j in range(n):
        minor = np.delete(np.delete(A, 0, axis=0), j, axis=1)
        cofactor = _tf_det(minor)
        term = _as_tf(A[0, j]) * cofactor
        out = out + term if (j % 2 == 0) else out - term
    return out


def poles_of_resolvent(M: np.ndarray, Delta: np.ndarray) -> np.ndarray:
    try:
        sys = feedback_resolvent(M, Delta)
        return np.asarray(control.poles(sys), dtype=complex)
    except ControlMIMONotImplemented:
        # Fallback when python-control lacks MIMO TF->SS support (e.g., no slycot):
        # poles of (I - MΔ)^(-1) are zeros of det(I - MΔ).
        L = tf_matmul(M, Delta)
        C = eye_tf(L.shape[0])
        for i in range(C.shape[0]):
            for j in range(C.shape[1]):
                C[i, j] = C[i, j] - _as_tf(L[i, j])
        det_char = _tf_det(C)
        return np.asarray(control.zeros(det_char), dtype=complex)


def is_unstable(M: np.ndarray, Delta: np.ndarray, pole_tol: float = 1e-8, near_tol: float = 1e-6) -> bool:
    """True if any pole has real part >= -pole_tol."""
    poles = poles_of_resolvent(M, Delta)
    if poles.size == 0:
        return False
    max_real = float(np.max(np.real(poles)))
    if abs(max_real) < near_tol:
        warnings.warn(
            f"Near-marginal closed-loop pole detected: max real part {max_real:.3e}",
            RuntimeWarning,
        )
    return bool(np.any(np.real(poles) >= -pole_tol))
