"""Publication-style plotting helpers for Ga₂O₃ defect analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MULLIKEN_SPIN_HSE: dict[str, dict[str, float]] = {
    "TETRA": {"O_I": 1.074, "O_II": 2.089, "O_III": 0.091},
    "OCTA": {"O_I": 2.009, "O_II": 0.172, "O_III": 1.053},
}
"""Spin populations in μB used for the oxygen-sublattice summary plot."""


def configure_matplotlib() -> None:
    """Apply a compact set of Matplotlib defaults used by the figures."""
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "figure.dpi": 150,
    })


def plot_dos(data: pd.DataFrame, *, energy_window: tuple[float, float] | None = None) -> tuple[plt.Figure, plt.Axes]:
    """Plot one or more DOS columns against energy in eV."""
    configure_matplotlib()
    y_columns = tuple(column for column in data.columns if column.startswith("y"))
    fig, ax = plt.subplots(figsize=(4.0, 4.2))
    for column in y_columns:
        ax.plot(data[column], data["energy_ev"], label=column)
    if energy_window is not None:
        ax.set_ylim(*energy_window)
    ax.axhline(0.0, linewidth=0.8, linestyle="--")
    ax.set_xlabel("DOS")
    ax.set_ylabel("Energy (eV)")
    if len(y_columns) > 1:
        ax.legend(frameon=False)
    return fig, ax


def plot_bands(bands: pd.DataFrame, *, energy_window: tuple[float, float] | None = None) -> tuple[plt.Figure, plt.Axes]:
    """Plot long-form band data with columns ``k``, ``band`` and ``energy_ev``."""
    configure_matplotlib()
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    for _, subset in bands.groupby("band"):
        ax.plot(subset["k"], subset["energy_ev"], linewidth=0.8)
    if energy_window is not None:
        ax.set_ylim(*energy_window)
    ax.axhline(0.0, linewidth=0.8, linestyle="--")
    ax.set_xlabel("k-path")
    ax.set_ylabel("Energy (eV)")
    return fig, ax


def plot_mulliken_summary(out: str | Path | None = None) -> tuple[plt.Figure, plt.Axes]:
    """Plot the oxygen-sublattice Mulliken spin population summary."""
    configure_matplotlib()
    labels = ("O_I", "O_II", "O_III")
    x = np.arange(len(labels), dtype=float)
    width = 0.34
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for offset, (case, values) in zip((-width / 2, width / 2), MULLIKEN_SPIN_HSE.items()):
        y = [values[label] for label in labels]
        ax.bar(x + offset, y, width=width, label=case, alpha=0.85)
        for xi, yi in zip(x + offset, y):
            ax.text(xi, yi + 0.04, f"{yi:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, [r"$O_{I}$", r"$O_{II}$", r"$O_{III}$"])
    ax.set_ylabel(r"Spin population ($\mu_B$)")
    ax.legend(frameon=False)
    if out is not None:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=600, bbox_inches="tight")
    return fig, ax
