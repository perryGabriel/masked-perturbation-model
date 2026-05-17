"""Bridge ANDES-derived state-space data into an MPM-compatible model map."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .channels import ChannelSelection
from .linearize import LinearModel


@dataclass(frozen=True)
class StateSpaceModelMap:
    """Minimal MPM-compatible model map for ``M(s)=C(sI-A)^{-1}B+D``.

    The standalone ANDES demo originally targeted a future/alternate
    ``mpmgame.lft.build_model_map`` API.  The current
    repository branch does not expose that import path, so this light wrapper
    keeps the verification notebook executable while preserving the same public
    dimensions and ``eval(s)`` behavior expected by the demo workflow.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray

    @property
    def nstates(self) -> int:
        return int(self.A.shape[0])

    @property
    def ninputs(self) -> int:
        return int(self.B.shape[1])

    @property
    def noutputs(self) -> int:
        return int(self.C.shape[0])

    def eval(self, s: complex) -> np.ndarray:
        """Evaluate the transfer matrix at the complex frequency ``s``."""
        identity = np.eye(self.nstates, dtype=complex)
        resolvent = np.linalg.solve(s * identity - self.A, self.B)
        return self.C @ resolvent + self.D


def _fallback_model_map(linear_model: LinearModel, channels: ChannelSelection) -> StateSpaceModelMap:
    return StateSpaceModelMap(
        A=np.asarray(linear_model.A, dtype=float),
        B=np.asarray(channels.B_w, dtype=float),
        C=np.asarray(channels.C_r, dtype=float),
        D=np.asarray(channels.D_rw, dtype=float),
    )


def build_lft_from_linear_model(linear_model: LinearModel, channels: ChannelSelection) -> Any:
    """Build an MPM-compatible model map ``M : w -> r`` from demo state-space data.

    If a future/full core API exposing
    ``mpmgame.lft.build_model_map`` is available, this helper
    will use it.  Otherwise it returns :class:`StateSpaceModelMap`, which
    supports the verification notebook's dimension checks and frequency-domain
    evaluation directly.
    """
    try:
        from mpmgame.lft import build_model_map
    except ImportError:
        return _fallback_model_map(linear_model, channels)

    return build_model_map(
        np.asarray(linear_model.A),
        np.asarray(channels.B_w),
        np.asarray(channels.C_r),
        np.asarray(channels.D_rw),
    )
