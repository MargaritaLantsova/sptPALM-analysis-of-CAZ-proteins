"""Plotting utilities for sptPALM analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def setup_axis(ax):
    """Apply common minimal plot styling."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.2)


def save_figure(
    fig,
    output_path: Path,
    dpi: int = 300,
    save_pdf: bool = False,
) -> None:
    """Save matplotlib figure and create parent directory if needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if save_pdf:
        fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")


def plot_track_length_distribution(
    spots: pd.DataFrame,
    min_track_length: int = 10,
    title: str = "Track length distribution",
    output_path: Path | None = None,
):
    """Plot track length distribution."""
    track_lengths = spots.groupby("track_id").size()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(track_lengths, bins=40)
    ax.axvline(
        min_track_length,
        linestyle="--",
        linewidth=2,
        label=f"min length = {min_track_length}",
    )

    ax.set_xlabel("Track length, frames")
    ax.set_ylabel("Number of tracks")
    ax.set_title(title)
    ax.legend(frameon=False)
    setup_axis(ax)
    plt.tight_layout()

    if output_path is not None:
        save_figure(fig, output_path)

    return fig, ax


def plot_edge_filter_example(
    spots: pd.DataFrame,
    bounds: pd.DataFrame,
    title: str = "Edge filter check",
    output_path: Path | None = None,
):
    """Visualize kept and removed localizations after edge filtering."""
    edge_ids = set(
        bounds.loc[bounds["is_edge_track"], "track_id"].to_numpy()
    )

    edge_spots = spots[spots["track_id"].isin(edge_ids)]
    kept_spots = spots[~spots["track_id"].isin(edge_ids)]

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        kept_spots["x"],
        kept_spots["y"],
        s=2,
        alpha=0.25,
        label="kept tracks",
    )

    ax.scatter(
        edge_spots["x"],
        edge_spots["y"],
        s=2,
        alpha=0.25,
        label="removed edge tracks",
    )

    ax.set_aspect("equal")
    ax.set_xlabel("x, µm")
    ax.set_ylabel("y, µm")
    ax.set_title(title)
    ax.legend(frameon=False, loc="center")
    setup_axis(ax)
    plt.tight_layout()

    if output_path is not None:
        save_figure(fig, output_path)

    return fig, ax


def plot_cell_paired_metric(
    data: pd.DataFrame,
    metric_column: str,
    protein: str,
    condition_order: list[str],
    condition_labels: dict[str, str],
    ylabel: str,
    title: str,
    output_path: Path | None = None,
):
    """Plot paired per-cell metric across conditions."""
    plot_data = data[data["protein"] == protein].copy()

    if plot_data.empty:
        return None, None

    plot_data["condition"] = pd.Categorical(
        plot_data["condition"],
        categories=condition_order,
        ordered=True,
    )

    plot_data = plot_data.sort_values(["sample", "cell", "condition"])

    x_map = {
        condition: index
        for index, condition in enumerate(condition_order)
    }

    fig, ax = plt.subplots(figsize=(7, 5))

    for _, group in plot_data.groupby(["sample", "cell"], observed=True):
        group = group.sort_values("condition")

        x = group["condition"].map(x_map).astype(float).to_numpy()
        y = group[metric_column].to_numpy()

        ax.plot(
            x,
            y,
            marker="o",
            linewidth=1.5,
            alpha=0.55,
        )

    median_summary = (
        plot_data.groupby("condition", observed=True)[metric_column]
        .median()
        .reindex(condition_order)
    )

    ax.plot(
        range(len(condition_order)),
        median_summary.values,
        marker="o",
        linewidth=4,
        color="black",
        label="group median",
    )

    ax.set_xticks(range(len(condition_order)))
    ax.set_xticklabels(
        [condition_labels[c] for c in condition_order],
        rotation=20,
    )

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    setup_axis(ax)
    ax.legend(frameon=False)
    plt.tight_layout()

    if output_path is not None:
        save_figure(fig, output_path)

    return fig, ax


def plot_ecdf(
    values_by_label: dict[str, np.ndarray],
    xlabel: str,
    title: str,
    xscale: str | None = None,
    output_path: Path | None = None,
):
    """Plot empirical cumulative distribution functions."""
    fig, ax = plt.subplots(figsize=(7, 5))

    for label, values in values_by_label.items():
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]

        if len(values) == 0:
            continue

        values = np.sort(values)
        cumulative = np.arange(1, len(values) + 1) / len(values)

        ax.plot(values, cumulative, linewidth=2, label=label)

    if xscale is not None:
        ax.set_xscale(xscale)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cumulative fraction")
    ax.set_title(title)
    ax.legend(frameon=False)
    setup_axis(ax)
    plt.tight_layout()

    if output_path is not None:
        save_figure(fig, output_path)

    return fig, ax
