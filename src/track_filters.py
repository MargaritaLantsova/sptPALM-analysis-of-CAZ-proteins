"""Trajectory filtering utilities for sptPALM analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd


def get_track_bounds(spots: pd.DataFrame) -> pd.DataFrame:
    """Calculate spatial bounds and length for each trajectory."""
    return (
        spots.groupby("track_id", sort=False)
        .agg(
            xmin=("x", "min"),
            xmax=("x", "max"),
            ymin=("y", "min"),
            ymax=("y", "max"),
            track_length=("frame", "count"),
            first_frame=("frame", "min"),
            last_frame=("frame", "max"),
        )
        .reset_index()
    )


def classify_edge_tracks(
    spots: pd.DataFrame,
    margin_um: float = 0.5,
) -> pd.DataFrame:
    """Classify tracks as edge or non-edge based on image-border margin."""
    x_min_image = spots["x"].min()
    x_max_image = spots["x"].max()
    y_min_image = spots["y"].min()
    y_max_image = spots["y"].max()

    bounds = get_track_bounds(spots)

    bounds["is_edge_track"] = (
        (bounds["xmin"] <= x_min_image + margin_um)
        | (bounds["xmax"] >= x_max_image - margin_um)
        | (bounds["ymin"] <= y_min_image + margin_um)
        | (bounds["ymax"] >= y_max_image - margin_um)
    )

    return bounds


def remove_edge_tracks(
    spots: pd.DataFrame,
    margin_um: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove tracks touching the image-border zone."""
    bounds = classify_edge_tracks(spots, margin_um=margin_um)

    edge_track_ids = bounds.loc[
        bounds["is_edge_track"],
        "track_id",
    ]

    filtered = spots[
        ~spots["track_id"].isin(edge_track_ids)
    ].copy()

    return filtered, bounds


def filter_short_tracks(
    spots: pd.DataFrame,
    min_track_length: int = 10,
) -> pd.DataFrame:
    """Remove tracks shorter than min_track_length."""
    track_lengths = spots.groupby("track_id").size()
    valid_track_ids = track_lengths[
        track_lengths >= min_track_length
    ].index

    return spots[spots["track_id"].isin(valid_track_ids)].copy()


def filter_tracks_by_displacement(
    spots: pd.DataFrame,
    min_total_displacement_um: float | None = None,
    max_total_displacement_um: float | None = None,
) -> pd.DataFrame:
    """Filter tracks by total displacement between first and last point."""
    rows = []

    for track_id, track in spots.groupby("track_id", sort=False):
        track = track.sort_values("frame")

        total_displacement = np.sqrt(
            (track["x"].iloc[-1] - track["x"].iloc[0]) ** 2
            + (track["y"].iloc[-1] - track["y"].iloc[0]) ** 2
        )

        rows.append(
            {
                "track_id": track_id,
                "total_displacement_um": total_displacement,
            }
        )

    displacement_table = pd.DataFrame(rows)

    mask = pd.Series(True, index=displacement_table.index)

    if min_total_displacement_um is not None:
        mask &= (
            displacement_table["total_displacement_um"]
            >= min_total_displacement_um
        )

    if max_total_displacement_um is not None:
        mask &= (
            displacement_table["total_displacement_um"]
            <= max_total_displacement_um
        )

    valid_ids = displacement_table.loc[mask, "track_id"]

    return spots[spots["track_id"].isin(valid_ids)].copy()


def preprocess_spots(
    spots: pd.DataFrame,
    edge_margin_um: float = 0.5,
    min_track_length: int = 10,
) -> tuple[pd.DataFrame, dict]:
    """Apply standard preprocessing: edge filtering and length filtering."""
    n_tracks_raw = spots["track_id"].nunique()
    n_spots_raw = len(spots)

    no_edge, edge_bounds = remove_edge_tracks(
        spots,
        margin_um=edge_margin_um,
    )

    filtered = filter_short_tracks(
        no_edge,
        min_track_length=min_track_length,
    )

    n_edge_tracks = int(edge_bounds["is_edge_track"].sum())

    qc = {
        "n_spots_raw": n_spots_raw,
        "n_tracks_raw": n_tracks_raw,
        "n_edge_tracks_removed": n_edge_tracks,
        "edge_tracks_removed_percent": (
            100 * n_edge_tracks / n_tracks_raw
            if n_tracks_raw
            else np.nan
        ),
        "n_tracks_after_edge": no_edge["track_id"].nunique(),
        "n_tracks_final": filtered["track_id"].nunique(),
        "edge_margin_um": edge_margin_um,
        "min_track_length": min_track_length,
    }

    return filtered, qc
