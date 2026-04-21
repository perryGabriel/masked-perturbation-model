"""Dynamic DSF parameter container for Q(z) transfer matrices.

This module is intentionally independent of optimizer internals. It provides a
small typed API to:
- hold DSF coefficient tensors,
- construct object-matrix transfer-function forms using :mod:`mpmgame.tf_tools`,
- enforce/validate hollow/mask/degree constraints,
- serialize/deserialize benchmark-row-safe payloads.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import control
import numpy as np

from .tf_tools import tf_matrix


def _as_bool_mask(mask: np.ndarray | None, shape: tuple[int, int]) -> np.ndarray:
    if mask is None:
        return np.ones(shape, dtype=bool)
    out = np.asarray(mask, dtype=bool)
    if out.shape != shape:
        raise ValueError(f"mask shape {out.shape} must equal matrix shape {shape}")
    return out


def _poly_desc_from_asc(coeffs: np.ndarray) -> np.ndarray:
    """Convert ascending-power coefficients to descending-power coefficients."""
    c = np.asarray(coeffs, dtype=float).reshape(-1)
    if c.size == 0:
        return np.array([0.0], dtype=float)
    return c[::-1]


@dataclass(slots=True)
class DynamicQ:
    """Typed container for dynamic Q(z) DSF entries.

    Coefficients are stored in ascending powers:
    - numerator: ``c0 + c1 z^-1 + ...``
    - denominator: ``1 + d1 z^-1 + ...`` (denominator coefficients exclude the
      leading 1)

    ``den_coeffs`` can be either:
    - shared denominator: shape ``(den_degree,)``
    - entry-specific denominator: shape ``(n, n, den_degree)``
    """

    num_coeffs: np.ndarray
    den_coeffs: np.ndarray | None = None
    mask: np.ndarray | None = None
    zero_diagonal: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_tensors(
        cls,
        num_coeffs: np.ndarray,
        den_coeffs: np.ndarray | None = None,
        mask: np.ndarray | None = None,
        *,
        zero_diagonal: bool = True,
        metadata: dict[str, Any] | None = None,
        enforce: bool = True,
    ) -> "DynamicQ":
        q = cls(
            num_coeffs=np.asarray(num_coeffs, dtype=float).copy(),
            den_coeffs=None if den_coeffs is None else np.asarray(den_coeffs, dtype=float).copy(),
            mask=None if mask is None else np.asarray(mask, dtype=bool).copy(),
            zero_diagonal=bool(zero_diagonal),
            metadata={} if metadata is None else dict(metadata),
        )
        q.validate()
        if enforce:
            q = q.enforced()
        return q

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(self.num_coeffs.shape[:2])

    @property
    def n(self) -> int:
        p, q = self.shape
        if p != q:
            raise ValueError(f"Q must be square; got {self.shape}")
        return p

    @property
    def num_degree(self) -> int:
        return int(self.num_coeffs.shape[2] - 1)

    @property
    def den_degree(self) -> int:
        if self.den_coeffs is None:
            return 0
        if self.den_coeffs.ndim == 1:
            return int(self.den_coeffs.shape[0])
        if self.den_coeffs.ndim == 3:
            return int(self.den_coeffs.shape[2])
        raise ValueError("den_coeffs must be 1D (shared) or 3D (entry-specific)")

    @property
    def uses_shared_denominator(self) -> bool:
        return self.den_coeffs is not None and self.den_coeffs.ndim == 1

    def validate(self) -> None:
        if self.num_coeffs.ndim != 3:
            raise ValueError("num_coeffs must have shape (n, n, num_degree+1)")
        n, m, terms = self.num_coeffs.shape
        if n != m:
            raise ValueError(f"Q must be square; got num_coeffs shape {self.num_coeffs.shape}")
        if terms < 1:
            raise ValueError("num_coeffs last axis must have at least one term")
        _ = _as_bool_mask(self.mask, (n, n))

        if self.den_coeffs is None:
            return
        if self.den_coeffs.ndim == 1:
            return
        if self.den_coeffs.ndim == 3 and self.den_coeffs.shape[:2] == (n, n):
            return
        raise ValueError(
            "den_coeffs must be None, shape (den_degree,), or shape (n, n, den_degree)"
        )

    def verify_constraints(self, tol: float = 1e-12) -> list[str]:
        """Return list of violated constraint tags."""
        self.validate()
        reasons: list[str] = []
        n = self.n
        mask = _as_bool_mask(self.mask, (n, n))

        if self.zero_diagonal:
            diag_num = self.num_coeffs[np.arange(n), np.arange(n), :]
            if np.any(np.abs(diag_num) > tol):
                reasons.append("non_hollow_numerator")
            if self.den_coeffs is not None and self.den_coeffs.ndim == 3:
                diag_den = self.den_coeffs[np.arange(n), np.arange(n), :]
                if np.any(np.abs(diag_den) > tol):
                    reasons.append("non_hollow_denominator")

        forbidden = (~mask) & (~np.eye(n, dtype=bool))
        if np.any(np.abs(self.num_coeffs[forbidden, :]) > tol):
            reasons.append("mask_violation_numerator")
        if self.den_coeffs is not None and self.den_coeffs.ndim == 3 and np.any(np.abs(self.den_coeffs[forbidden, :]) > tol):
            reasons.append("mask_violation_denominator")

        if self.num_degree < 0 or self.den_degree < 0:
            reasons.append("negative_degree")

        return reasons

    def enforced(self) -> "DynamicQ":
        """Return a copy with diagonal/forbidden entries zeroed by constraints."""
        self.validate()
        n = self.n
        mask = _as_bool_mask(self.mask, (n, n))
        num = self.num_coeffs.copy()
        den = None if self.den_coeffs is None else self.den_coeffs.copy()

        forbidden = (~mask) & (~np.eye(n, dtype=bool))
        num[forbidden, :] = 0.0
        if den is not None and den.ndim == 3:
            den[forbidden, :] = 0.0

        if self.zero_diagonal:
            num[np.arange(n), np.arange(n), :] = 0.0
            if den is not None and den.ndim == 3:
                den[np.arange(n), np.arange(n), :] = 0.0

        return DynamicQ(
            num_coeffs=num,
            den_coeffs=den,
            mask=mask,
            zero_diagonal=self.zero_diagonal,
            metadata=dict(self.metadata),
        )

    def to_tf_matrix(self) -> np.ndarray:
        """Build an object matrix of SISO transfer functions."""
        q = self.enforced()
        n = q.n
        entries: list[list[control.TransferFunction]] = []
        for i in range(n):
            row: list[control.TransferFunction] = []
            for j in range(n):
                num_desc = _poly_desc_from_asc(q.num_coeffs[i, j, :])
                if q.den_coeffs is None:
                    den_asc = np.array([1.0], dtype=float)
                elif q.den_coeffs.ndim == 1:
                    den_asc = np.concatenate([[1.0], q.den_coeffs])
                else:
                    den_asc = np.concatenate([[1.0], q.den_coeffs[i, j, :]])
                den_desc = _poly_desc_from_asc(den_asc)
                row.append(control.tf(num_desc, den_desc))
            entries.append(row)
        return tf_matrix(entries)

    def to_benchmark_payload(self) -> dict[str, Any]:
        """JSON-safe structure suitable for benchmark-row serialization."""
        payload: dict[str, Any] = {
            "kind": "dsf_dynamic_q",
            "shape": list(self.shape),
            "num_degree": self.num_degree,
            "den_degree": self.den_degree,
            "shared_denominator": self.uses_shared_denominator,
            "zero_diagonal": self.zero_diagonal,
            "num_coeffs": self.num_coeffs.tolist(),
            "den_coeffs": None if self.den_coeffs is None else self.den_coeffs.tolist(),
            "mask": None if self.mask is None else np.asarray(self.mask, dtype=bool).astype(int).tolist(),
            "metadata": dict(self.metadata),
        }
        return payload

    @classmethod
    def from_benchmark_payload(cls, payload: dict[str, Any]) -> "DynamicQ":
        kind = str(payload.get("kind", ""))
        if kind not in {"", "dsf_dynamic_q"}:
            raise ValueError(f"Unsupported payload kind={kind}")
        mask_obj = payload.get("mask", None)
        return cls.from_tensors(
            num_coeffs=np.asarray(payload["num_coeffs"], dtype=float),
            den_coeffs=None if payload.get("den_coeffs", None) is None else np.asarray(payload["den_coeffs"], dtype=float),
            mask=None if mask_obj is None else np.asarray(mask_obj, dtype=bool),
            zero_diagonal=bool(payload.get("zero_diagonal", True)),
            metadata=dict(payload.get("metadata", {})),
            enforce=False,
        )


def serialize_dynamic_q_for_row(dynamic_q: DynamicQ) -> dict[str, Any]:
    """Helper alias for benchmark codepaths that build per-row payloads."""
    return dynamic_q.to_benchmark_payload()


def deserialize_dynamic_q_from_row(payload: dict[str, Any]) -> DynamicQ:
    """Inverse helper for benchmark-row payloads."""
    return DynamicQ.from_benchmark_payload(payload)
