"""Spatial clustering utilities for localization data."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


def run_dbscan_on_spots(
    spots: pd.DataFrame,
    eps_um: float = 0.3,
    min_samples: int = 10,
) -> pd.DataFrame:
    """Run DBSCAN clustering on x/y localization coordinates."""
    if spots.empty:
        result = spots.copy()
        result["cluster"] = []
        return result

    coords = spots[["x", "y"]].to_numpy(dtype=float)

    clustering = DBSCAN(
        eps=eps_um,
        min_samples=min_samples,
    ).fit(coords)

    result = spots.copy()
    result["cluster"] = clustering.labels_

    return result


def summarize_dbscan_clusters(clustered_spots: pd.DataFrame) -> dict:
    """Summarize DBSCAN clustering result."""
    if clustered_spots.empty or "cluster" not in clustered_spots.columns:
        return {
            "n_clusters": 0,
            "clustered_fraction": 0,
            "median_cluster_size": np.nan,
            "mean_cluster_size": np.nan,
            "n_noise_points": 0,
        }

    clustered = clustered_spots[clustered_spots["cluster"] != -1].copy()
    noise = clustered_spots[clustered_spots["cluster"] == -1].copy()

    n_total = len(clustered_spots)
    n_clustered = len(clustered)

    if n_clustered == 0:
        return {
            "n_clusters": 0,
            "clustered_fraction": 0,
            "median_cluster_size": np.nan,
            "mean_cluster_size": np.nan,
            "n_noise_points": len(noise),
        }

    cluster_sizes = clustered.groupby("cluster").size()

    return {
        "n_clusters": cluster_sizes.shape[0],
        "clustered_fraction": n_clustered / n_total,
        "median_cluster_size": cluster_sizes.median(),
        "mean_cluster_size": cluster_sizes.mean(),
        "n_noise_points": len(noise),
    }


def batch_dbscan_summary(
    spots: pd.DataFrame,
    eps_um: float = 0.3,
    min_samples: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run DBSCAN per file/cell and return clustered spots and summary."""
    clustered_tables = []
    summary_rows = []

    group_columns = [
        "protein",
        "condition",
        "sample",
        "cell",
        "file_name",
    ]

    for keys, group in spots.groupby(group_columns, sort=False):
        protein, condition, sample, cell, file_name = keys

        clustered = run_dbscan_on_spots(
            group,
            eps_um=eps_um,
            min_samples=min_samples,
        )

        summary = summarize_dbscan_clusters(clustered)

        summary_rows.append(
            {
                "protein": protein,
                "condition": condition,
                "sample": sample,
                "cell": cell,
                "file_name": file_name,
                "eps_um": eps_um,
                "min_samples": min_samples,
                **summary,
            }
        )

        clustered_tables.append(clustered)

    clustered_spots = (
        pd.concat(clustered_tables, ignore_index=True)
        if clustered_tables
        else pd.DataFrame()
    )

    summary_table = pd.DataFrame(summary_rows)

    return clustered_spots, summary_table


def calculate_cluster_centroids(clustered_spots: pd.DataFrame) -> pd.DataFrame:
    """Calculate centroid and size for each DBSCAN cluster."""
    if clustered_spots.empty or "cluster" not in clustered_spots.columns:
        return pd.DataFrame()

    clustered = clustered_spots[clustered_spots["cluster"] != -1].copy()

    if clustered.empty:
        return pd.DataFrame()

    group_columns = [
        "protein",
        "condition",
        "sample",
        "cell",
        "file_name",
        "cluster",
    ]

    return (
        clustered.groupby(group_columns, observed=True)
        .agg(
            centroid_x_um=("x", "mean"),
            centroid_y_um=("y", "mean"),
            cluster_size=("cluster", "count"),
        )
        .reset_index()
    )
