"""Edge-artifact analysis utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .track_filters import classify_edge_tracks, remove_edge_tracks


def calculate_edge_summary(
    spots: pd.DataFrame,
    margin_um: float = 0.5,
) -> dict:
    """Calculate edge-filtering summary for one field of view."""
    bounds = classify_edge_tracks(spots, margin_um=margin_um)

    n_tracks_raw = bounds["track_id"].nunique()
    n_edge_tracks = int(bounds["is_edge_track"].sum())

    return {
        "n_tracks_raw": n_tracks_raw,
        "n_edge_tracks_removed": n_edge_tracks,
        "edge_tracks_removed_percent": (
            100 * n_edge_tracks / n_tracks_raw if n_tracks_raw else np.nan
        ),
        "edge_margin_um": margin_um,
    }


def edge_sensitivity_for_spots(
    spots: pd.DataFrame,
    margins_um: list[float] | tuple[float, ...],
) -> pd.DataFrame:
    """Calculate edge removal fraction across multiple margins."""
    rows = []

    for margin in margins_um:
        summary = calculate_edge_summary(spots, margin_um=margin)
        rows.append(summary)

    return pd.DataFrame(rows)


def compare_edge_vs_center_links(
    spots: pd.DataFrame,
    margin_um: float = 0.5,
) -> pd.DataFrame:
    """Compare displacement values for links near edge and image center."""
    x_min = spots["x"].min()
    x_max = spots["x"].max()
    y_min = spots["y"].min()
    y_max = spots["y"].max()

    data = spots.copy()
    data["is_edge_point"] = (
        (data["x"] <= x_min + margin_um)
        | (data["x"] >= x_max - margin_um)
        | (data["y"] <= y_min + margin_um)
        | (data["y"] >= y_max - margin_um)
    )

    rows = []

    for track_id, track in data.groupby("track_id", sort=False):
        track = track.sort_values("frame")

        if len(track) < 2:
            continue

        x = track["x"].to_numpy(dtype=float)
        y = track["y"].to_numpy(dtype=float)
        edge = track["is_edge_point"].to_numpy(dtype=bool)

        displacement = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
        link_is_edge = edge[:-1] | edge[1:]

        for value, is_edge in zip(displacement, link_is_edge):
            rows.append(
                {
                    "track_id": track_id,
                    "displacement_um": value,
                    "region": "edge" if is_edge else "center",
                }
            )

    return pd.DataFrame(rows)


def calculate_edge_filtered_tables(
    spots: pd.DataFrame,
    margin_um: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return raw spots, edge-filtered spots, and edge bounds table."""
    filtered_spots, bounds = remove_edge_tracks(
        spots,
        margin_um=margin_um,
    )

    return spots.copy(), filtered_spots, bounds
