"""Diffusion analysis functions for single-particle trajectories."""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_msd_for_track(track: pd.DataFrame) -> pd.DataFrame:
    """Compute time-averaged MSD for one trajectory."""
    track = track.sort_values("frame")

    x = track["x"].to_numpy(dtype=float)
    y = track["y"].to_numpy(dtype=float)
    t = track["t"].to_numpy(dtype=float)
    frame = track["frame"].to_numpy(dtype=int)

    n_points = len(track)

    if n_points < 2:
        return pd.DataFrame()

    rows = []

    for lag in range(1, n_points):
        dx = x[lag:] - x[:-lag]
        dy = y[lag:] - y[:-lag]

        squared_displacement = dx**2 + dy**2

        dt = t[lag:] - t[:-lag]
        tau_s = np.nanmedian(dt)

        if not np.isfinite(tau_s) or tau_s <= 0:
            frame_dt = frame[lag:] - frame[:-lag]
            tau_s = np.nanmedian(frame_dt)

        rows.append(
            {
                "lag_frames": lag,
                "tau_s": tau_s,
                "msd_um2": np.nanmean(squared_displacement),
                "n_displacements": len(squared_displacement),
            }
        )

    return pd.DataFrame(rows)


def estimate_diffusion_coefficient(
    msd: pd.DataFrame,
    n_fit_points: int = 3,
    min_msd_points: int = 3,
) -> tuple[float, float, float, float]:
    """Estimate diffusion coefficient from initial MSD points.

    For 2D Brownian diffusion, MSD(t) = 4Dt, so D = slope / 4.
    """
    if msd.empty or len(msd) < min_msd_points:
        return np.nan, np.nan, np.nan, np.nan

    fit_data = msd.dropna(subset=["tau_s", "msd_um2"]).copy()
    fit_data = fit_data[
        (fit_data["tau_s"] > 0) & (fit_data["msd_um2"] >= 0)
    ].head(n_fit_points)

    if len(fit_data) < min_msd_points:
        return np.nan, np.nan, np.nan, np.nan

    tau = fit_data["tau_s"].to_numpy(dtype=float)
    msd_values = fit_data["msd_um2"].to_numpy(dtype=float)

    slope, intercept = np.polyfit(tau, msd_values, deg=1)
    diffusion_coefficient = slope / 4

    predicted = slope * tau + intercept
    residual_sum_of_squares = np.sum((msd_values - predicted) ** 2)
    total_sum_of_squares = np.sum((msd_values - np.mean(msd_values)) ** 2)

    r_squared = (
        1 - residual_sum_of_squares / total_sum_of_squares
        if total_sum_of_squares > 0
        else np.nan
    )

    return diffusion_coefficient, slope, intercept, r_squared


def calculate_track_diffusion(
    spots: pd.DataFrame,
    n_fit_points: int = 3,
    min_msd_points: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate MSD and diffusion coefficient for all trajectories."""
    track_results = []
    msd_results = []

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

        msd = compute_msd_for_track(track)

        if msd.empty:
            continue

        diffusion_coefficient, slope, intercept, r_squared = (
            estimate_diffusion_coefficient(
                msd,
                n_fit_points=n_fit_points,
                min_msd_points=min_msd_points,
            )
        )

        track = track.sort_values("frame")

        total_displacement = np.sqrt(
            (track["x"].iloc[-1] - track["x"].iloc[0]) ** 2
            + (track["y"].iloc[-1] - track["y"].iloc[0]) ** 2
        )

        track_results.append(
            {
                "protein": protein,
                "condition": condition,
                "sample": sample,
                "cell": cell,
                "file_name": file_name,
                "track_id": track_id,
                "track_length": len(track),
                "total_displacement_um": total_displacement,
                "D_um2_s": diffusion_coefficient,
                "log10_D": (
                    np.log10(diffusion_coefficient)
                    if np.isfinite(diffusion_coefficient)
                    and diffusion_coefficient > 0
                    else np.nan
                ),
                "msd_slope": slope,
                "msd_intercept": intercept,
                "msd_fit_r2": r_squared,
                "n_fit_points": n_fit_points,
            }
        )

        msd = msd.copy()
        msd["protein"] = protein
        msd["condition"] = condition
        msd["sample"] = sample
        msd["cell"] = cell
        msd["file_name"] = file_name
        msd["track_id"] = track_id

        msd_results.append(msd)

    track_results = pd.DataFrame(track_results)
    msd_results = (
        pd.concat(msd_results, ignore_index=True)
        if msd_results
        else pd.DataFrame()
    )

    return track_results, msd_results


def summarize_diffusion_by_cell(track_results: pd.DataFrame) -> pd.DataFrame:
    """Create per-cell diffusion summary table."""
    valid = track_results[
        np.isfinite(track_results["D_um2_s"])
        & (track_results["D_um2_s"] > 0)
    ].copy()

    if valid.empty:
        return pd.DataFrame()

    return (
        valid.groupby(
            ["protein", "condition", "sample", "cell", "file_name"],
            observed=True,
        )
        .agg(
            n_tracks=("track_id", "count"),
            median_D_um2_s=("D_um2_s", "median"),
            mean_D_um2_s=("D_um2_s", "mean"),
            median_log10_D=("log10_D", "median"),
            mean_log10_D=("log10_D", "mean"),
            median_track_length=("track_length", "median"),
            median_total_displacement_um=("total_displacement_um", "median"),
        )
        .reset_index()
    )
