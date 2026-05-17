"""Bridge ANDES-derived state-space data into the core MPM LFT utilities."""
from __future__ import annotations

import numpy as np

from mpmgame.lft import build_model_map

from .channels import ChannelSelection
from .linearize import LinearModel


def build_lft_from_linear_model(linear_model: LinearModel, channels: ChannelSelection):
    """Build the LFT-derived MPM model map ``M : w -> r``.

    The ANDES/grid demo is responsible for producing the channel realization

    ``x_dot = A x + B_w w``
    ``r = C_r x + D_rw w``.

    The reusable MPM package owns the LFT/model-map construction itself, via
    :func:`mpmgame.lft.build_model_map`. The returned object is the transfer-
    function matrix ``M(s)=C_r(sI-A)^(-1)B_w + D_rw`` used directly by the
    existing MPM attack/defense/game routines.
    """
    return build_model_map(
        np.asarray(linear_model.A, dtype=float),
        np.asarray(channels.B_w, dtype=float),
        np.asarray(channels.C_r, dtype=float),
        np.asarray(channels.D_rw, dtype=float),
    )
