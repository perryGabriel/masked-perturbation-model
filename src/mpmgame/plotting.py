"""Optional visualization utilities for MPM games."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_payoff_heatmap(payoff: np.ndarray, attack_labels: list[str], defense_labels: list[str]):
    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(payoff, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(defense_labels)), defense_labels)
    ax.set_yticks(range(len(attack_labels)), attack_labels)
    ax.set_xlabel("Defense")
    ax.set_ylabel("Attack")
    ax.set_title("Attacker payoff")
    fig.colorbar(im, ax=ax)
    return fig, ax


def plot_success_bipartite_graph(success_sets: dict[str, set[str]]):
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError("networkx is required for plot_success_bipartite_graph") from exc

    G = nx.Graph()
    attacks = list(success_sets.keys())
    defenses = sorted(set().union(*success_sets.values())) if success_sets else []
    G.add_nodes_from(attacks, bipartite=0)
    G.add_nodes_from(defenses, bipartite=1)
    for a, ds in success_sets.items():
        for d in ds:
            G.add_edge(a, d)

    pos = {}
    for i, a in enumerate(attacks):
        pos[a] = (0, -i)
    for j, d in enumerate(defenses):
        pos[d] = (1, -j)

    fig, ax = plt.subplots(figsize=(6, 4))
    nx.draw(G, pos=pos, with_labels=True, ax=ax, node_size=1200)
    ax.set_title("Attack-success bipartite graph")
    return fig, ax
