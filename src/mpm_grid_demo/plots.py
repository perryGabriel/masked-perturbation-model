"""Plotting helpers for the ANDES verification demo."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_eigs(eigvals: np.ndarray, title: str, path: str | Path | None = None):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(np.real(eigvals), np.imag(eigvals))
    ax.axvline(0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("Real part")
    ax.set_ylabel("Imaginary part")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if path is not None:
        fig.savefig(path, dpi=200)
    return fig, ax


def plot_heatmap(
    V: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    path: str | Path | None = None,
):
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(V, aspect="auto")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    if path is not None:
        fig.savefig(path, dpi=200)
    return fig, ax


def plot_game_matrix(df: pd.DataFrame, path: str | Path | None = None):
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * df.shape[1]), max(4, 0.3 * df.shape[0])))
    im = ax.imshow(df.values, aspect="auto", vmin=0, vmax=1)
    ax.set_title("Attack-defense utility matrix")
    ax.set_xlabel("Defense mask")
    ax.set_ylabel("Attack")
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns, rotation=90)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    if path is not None:
        fig.savefig(path, dpi=200)
    return fig, ax
