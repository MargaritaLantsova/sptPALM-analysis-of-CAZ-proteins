"""Velocity and frame-to-frame displacement analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd


def calculate_frame_to_frame_links(spots: pd.DataFrame) -> pd.DataFrame:
    """Calculate frame-to-frame displacement and velocity for all trajectories."""
    link_rows = []

    group_columns = [
        "protein",
        "condition",
        "sample",
        "cell",
        "file_name",
        "track_id",
    ]

    for keys, track in spots.groupby(group_columns, sort=False):
        protein, condition, sample, cell, file_name, track_id = keys
        track = track.sort_values("frame")

        if len(track) < 2:
            continue

        x = track["x"].to_numpy(dtype=float)
        y = track["y"].to_numpy(dtype=float)
        t = track["t"].to_numpy(dtype=float)
        frames = track["frame"].to_numpy(dtype=int)

        dx = np.diff(x)
        dy = np.diff(y)
        dt = np.diff(t)
        frame_step = np.diff(frames)

        displacement_um = np.sqrt(dx**2 + dy**2)

        with np.errstate(divide="ignore", invalid="ignore"):
            velocity_um_s = displacement_um / dt

        for index in range(len(displacement_um)):
            if (
                not np.isfinite(velocity_um_s[index])
                or velocity_um_s[index] <= 0
            ):
                continue

            link_rows.append(
                {
                    "protein": protein,
                    "condition": condition,
                    "sample": sample,
                    "cell": cell,
                    "file_name": file_name,
                    "track_id": track_id,
                    "start_frame": frames[index],
                    "end_frame": frames[index + 1],
                    "frame_step": frame_step[index],
                    "dt_s": dt[index],
                    "dx_um": dx[index],
                    "dy_um": dy[index],
                    "displacement_um": displacement_um[index],
                    "velocity_um_s": velocity_um_s[index],
                }
            )

    return pd.DataFrame(link_rows)


def estimate_rayleigh_sigma(values: np.ndarray) -> float:
    """Estimate Rayleigh scale parameter from positive values."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]

    if len(values) == 0:
        return np.nan

    return np.sqrt(np.mean(values**2) / 2)


def rayleigh_pdf(x: np.ndarray, sigma: float) -> np.ndarray:
    """Rayleigh probability density function without scipy dependency."""
    x = np.asarray(x, dtype=float)

    if not np.isfinite(sigma) or sigma <= 0:
        return np.full_like(x, np.nan)

    return (x / sigma**2) * np.exp(-(x**2) / (2 * sigma**2))


def rayleigh_percentile(sigma: float, percentile: float) -> float:
    """Rayleigh percentile without scipy dependency."""
    probability = percentile / 100

    if not np.isfinite(sigma) or sigma <= 0:
        return np.nan

    return sigma * np.sqrt(-2 * np.log(1 - probability))


def calculate_control_velocity_thresholds(
    links: pd.DataFrame,
    percentile: float = 95,
) -> pd.DataFrame:
    """Calculate empirical and Rayleigh velocity thresholds from controls."""
    rows = []

    for protein in sorted(links["protein"].dropna().unique()):
        control_values = links.loc[
            (links["protein"] == protein)
            & (links["condition"] == "control"),
            "velocity_um_s",
        ].dropna().to_numpy()

        control_values = control_values[
            np.isfinite(control_values) & (control_values > 0)
        ]

        if len(control_values) == 0:
            continue

        empirical_threshold = np.percentile(control_values, percentile)
        sigma = estimate_rayleigh_sigma(control_values)
        rayleigh_threshold = rayleigh_percentile(sigma, percentile)

        rows.append(
            {
                "protein": protein,
                "threshold_percentile": percentile,
                "empirical_threshold_um_s": empirical_threshold,
                "rayleigh_sigma_um_s": sigma,
                "rayleigh_threshold_um_s": rayleigh_threshold,
                "n_control_links": len(control_values),
            }
        )

    return pd.DataFrame(rows)


def summarize_velocity_by_cell(
    links: pd.DataFrame,
    velocity_thresholds: pd.DataFrame,
) -> pd.DataFrame:
    """Create per-cell velocity summary using control-derived thresholds."""
    rows = []

    for keys, group in links.groupby(
        ["protein", "condition", "sample", "cell", "file_name"],
        observed=True,
    ):
        protein, condition, sample, cell, file_name = keys

        threshold_row = velocity_thresholds[
            velocity_thresholds["protein"] == protein
        ]

        if threshold_row.empty:
            empirical_threshold = np.nan
            rayleigh_threshold = np.nan
        else:
            empirical_threshold = threshold_row[
                "empirical_threshold_um_s"
            ].iloc[0]
            rayleigh_threshold = threshold_row[
                "rayleigh_threshold_um_s"
            ].iloc[0]

        velocities = group["velocity_um_s"].dropna().to_numpy()
        displacements = group["displacement_um"].dropna().to_numpy()

        rows.append(
            {
                "protein": protein,
                "condition": condition,
                "sample": sample,
                "cell": cell,
                "file_name": file_name,
                "n_links": len(group),
                "median_displacement_um": np.median(displacements),
                "mean_displacement_um": np.mean(displacements),
                "median_velocity_um_s": np.median(velocities),
                "mean_velocity_um_s": np.mean(velocities),
                "velocity_threshold_um_s": empirical_threshold,
                "rayleigh_threshold_um_s": rayleigh_threshold,
                "high_velocity_fraction": (
                    np.mean(velocities > empirical_threshold)
                    if np.isfinite(empirical_threshold)
                    else np.nan
                ),
                "high_velocity_percent": (
                    100 * np.mean(velocities > empirical_threshold)
                    if np.isfinite(empirical_threshold)
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(rows)
