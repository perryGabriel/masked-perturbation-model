"""Read/write channel selection for ANDES-derived linear models."""
from __future__ import annotations

from dataclasses import dataclass
import random

import numpy as np


@dataclass
class ChannelSelection:
    B_w: np.ndarray
    C_r: np.ndarray
    D_rw: np.ndarray
    read_metadata: list[dict]
    write_metadata: list[dict]


def _find_speed_indices(state_names: list[str]) -> list[int]:
    tokens = ("omega", "speed", "w_")
    return [k for k, name in enumerate(state_names) if any(t in name.lower() for t in tokens)]


def _find_power_indices(state_names: list[str]) -> list[int]:
    tokens = ("pm", "pe", "power")
    return [k for k, name in enumerate(state_names) if any(t in name.lower() for t in tokens)]


def select_speed_reads_and_power_writes(A: np.ndarray, state_names: list[str]) -> ChannelSelection:
    """Expose speed-like reads and power-like additive write channels.

    If inertia-aware write gains are not available, the demo uses unit-gain
    injections into identified power-like state coordinates. If no power-like
    states are found, the returned write channel matrix has zero columns; the
    notebook then switches to ``select_random_channels`` as a smoke-test path.
    """
    n_x = A.shape[0]
    speed_indices = _find_speed_indices(state_names)
    power_indices = _find_power_indices(state_names)

    if not speed_indices:
        speed_indices = list(range(n_x // 2, n_x))

    n_r = len(speed_indices)
    n_w = len(power_indices)

    C_r = np.zeros((n_r, n_x))
    for row, idx in enumerate(speed_indices):
        C_r[row, idx] = 1.0

    B_w = np.zeros((n_x, n_w))
    for col, idx in enumerate(power_indices):
        B_w[idx, col] = 1.0

    D_rw = np.zeros((n_r, n_w))

    read_metadata = [
        {"index": row, "state_index": idx, "name": state_names[idx], "type": "speed_read"}
        for row, idx in enumerate(speed_indices)
    ]
    write_metadata = [
        {
            "index": col,
            "state_index": idx,
            "name": f"Pm_attack_to_{state_names[idx]}",
            "type": "power_write",
        }
        for col, idx in enumerate(power_indices)
    ]

    return ChannelSelection(
        B_w=B_w,
        C_r=C_r,
        D_rw=D_rw,
        read_metadata=read_metadata,
        write_metadata=write_metadata,
    )


def select_random_channels(
    A: np.ndarray,
    state_names: list[str],
    n_channels: int,
    random_seed: int,
) -> ChannelSelection:
    """Select disjoint random read and write state coordinates for smoke tests."""
    random.seed(random_seed)
    n_x = A.shape[0]
    if n_channels < 1:
        raise ValueError("n_channels must be positive")
    if 2 * n_channels > n_x:
        raise ValueError("2 * n_channels cannot exceed the number of states")

    selected_indices = random.sample(range(n_x), 2 * n_channels)
    selected_r_indices = selected_indices[:n_channels]
    selected_w_indices = selected_indices[n_channels:]

    C_r = np.zeros((n_channels, n_x))
    for row, idx in enumerate(selected_r_indices):
        C_r[row, idx] = 1.0

    B_w = np.zeros((n_x, n_channels))
    for col, idx in enumerate(selected_w_indices):
        B_w[idx, col] = 1.0

    D_rw = np.zeros((n_channels, n_channels))

    read_metadata = [
        {"index": row, "state_index": idx, "name": state_names[idx], "type": "random_read"}
        for row, idx in enumerate(selected_r_indices)
    ]
    write_metadata = [
        {
            "index": col,
            "state_index": idx,
            "name": f"random_attack_to_{state_names[idx]}",
            "type": "random_write",
        }
        for col, idx in enumerate(selected_w_indices)
    ]

    return ChannelSelection(
        B_w=B_w,
        C_r=C_r,
        D_rw=D_rw,
        read_metadata=read_metadata,
        write_metadata=write_metadata,
    )
