"""Bridge ANDES-derived state-space data into the core MPM package."""
from __future__ import annotations

from typing import Any

import numpy as np

from .channels import ChannelSelection
from .linearize import LinearModel


def build_lft_from_linear_model(linear_model: LinearModel, channels: ChannelSelection) -> Any:
    """Build the core MPM model map ``M : w -> r`` from demo state-space data.

    The import is intentionally lazy so the demo package remains importable in
    environments that have not installed the full ``masked_perturbation_model``
    API surface yet.
    """
    try:
        from masked_perturbation_model.lft import build_model_map
    except ImportError as exc:
        raise ImportError(
            "build_lft_from_linear_model requires "
            "masked_perturbation_model.lft.build_model_map"
        ) from exc

    return build_model_map(
        np.asarray(linear_model.A),
        np.asarray(channels.B_w),
        np.asarray(channels.C_r),
        np.asarray(channels.D_rw),
    )
