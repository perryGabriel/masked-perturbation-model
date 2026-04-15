"""Plotting helpers for chapter-2 experiments."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_vulnerability_vs_parameter(df: pd.DataFrame, x: str, title: str, threat_model: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df[x], df["measured"], label=f"Measured V ({threat_model})", linewidth=2)
    for col, label in [("lower_single", "Lower bound (single-link assumptions)"), ("lower_full", "Lower bound (full assumptions)"), ("upper_full", "Upper bound (full, Q=0 helper)")]:
        if col in df.columns:
            ax.plot(df[x], df[col], "--", label=label)
    ax.set_xlabel(x)
    ax.set_ylabel("vulnerability / bound")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig, ax


def plot_bound_scatter(df: pd.DataFrame, bound_col: str = "lower_bound", title: str = "Measured vs bound"):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(df[bound_col], df["measured"], alpha=0.7)
    lo = min(df[bound_col].min(), df["measured"].min())
    hi = max(df[bound_col].max(), df["measured"].max())
    ax.plot([lo, hi], [lo, hi], "k--", label="y=x")
    ax.set_xlabel(bound_col)
    ax.set_ylabel("measured")
    ax.set_title(title)
    ax.legend()
    return fig, ax


def plot_projection_convergence(iterations_df: pd.DataFrame):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(iterations_df["iter"], iterations_df["surrogate_obj"], marker="o")
    ax[0].set_title("Projection surrogate objective")
    ax[0].set_xlabel("iteration")
    ax[0].set_ylabel("objective")
    ax[1].plot(iterations_df["iter"], iterations_df["measured_vulnerability"], marker="o")
    ax[1].set_title("Measured vulnerability")
    ax[1].set_xlabel("iteration")
    ax[1].set_ylabel("V")
    for a in ax:
        a.grid(True, alpha=0.3)
    return fig, ax


def plot_design_comparison(labels: list[str], values: list[float], title: str = "Design comparison"):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, values)
    ax.set_ylabel("true measured vulnerability")
    ax.set_title(title)
    return fig, ax


def plot_access_model_comparison(df: pd.DataFrame, threat_model: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["access_model"], df["measured"], label=f"Measured V ({threat_model})")
    if "lower_full" in df.columns:
        ax.plot(df["access_model"], df["lower_full"], "ro--", label="Lower full bound")
    ax.set_ylabel("vulnerability / bound")
    ax.set_title(f"Access model comparison ({threat_model})")
    ax.legend()
    return fig, ax
