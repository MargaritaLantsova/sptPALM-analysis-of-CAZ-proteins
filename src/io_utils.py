"""Input/output utilities for TrackMate-based sptPALM analysis."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd


COLUMN_ALIASES = {
    "TRACK_ID": "track_id",
    "POSITION_X": "x",
    "POSITION_Y": "y",
    "POSITION_T": "t",
    "FRAME": "frame",
    "QUALITY": "quality",
    "SNR_CH1": "snr",
    "SIGNAL_NOISE_RATIO_CH1": "snr",
    "MEAN_INTENSITY_CH1": "mean_intensity",
}


def ensure_directories(paths: Iterable[Path]) -> None:
    """Create directories if they do not exist."""
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)


def standardize_column_name(column_name: str) -> str:
    """Convert TrackMate column names to a compact uppercase format."""
    column_name = str(column_name).strip()
    column_name = column_name.replace(" ", "_")
    column_name = column_name.replace("-", "_")
    column_name = column_name.replace("/", "_")
    column_name = re.sub(r"[^A-Za-z0-9_]+", "", column_name)
    column_name = re.sub(r"_+", "_", column_name)
    return column_name.upper()


def normalize_condition_name(text: str) -> str:
    """Normalize condition names from folders or file names."""
    text = str(text).lower()
    text = text.replace(" ", "_").replace("-", "_")

    if "30" in text and "h2o2" in text:
        return "500uM_H2O2_30min"
    if "500" in text and ("h2o2" in text or "ox" in text):
        return "500uM_H2O2"
    if "control" in text or "ctrl" in text:
        return "control"

    return "unknown"


def infer_protein(path: Path) -> str:
    """Infer protein name from path parts or file name."""
    full_text = " ".join(path.parts).lower()

    if "snap" in full_text:
        return "SNAP25"
    if "syx" in full_text or "syntaxin" in full_text:
        return "Syx"

    return "unknown"


def infer_condition(path: Path) -> str:
    """Infer condition from parent folders and file name."""
    for part in reversed(path.parts):
        condition = normalize_condition_name(part)
        if condition != "unknown":
            return condition
    return "unknown"


def infer_sample(path: Path) -> float | int:
    """Infer sample number from path."""
    full_text = " ".join(path.parts)

    patterns = [
        r"Sample[_\s-]*(\d+)",
        r"sample[_\s-]*(\d+)",
        r"C(\d+)[-_]?cs",
        r"cs(\d+)",
    ]

    for pattern in patterns:
        match = re.search(pattern, full_text)
        if match:
            return int(match.group(1))

    return np.nan


def infer_cell(path: Path) -> float | int:
    """Infer cell number from file name."""
    match = re.search(r"cell[_\s-]?(\d+)", path.name, flags=re.IGNORECASE)
    if match:
        return int(match.group(1))
    return np.nan


def parse_file_metadata(path: Path) -> dict:
    """Extract metadata from a TrackMate CSV file path."""
    path = Path(path)

    return {
        "protein": infer_protein(path),
        "condition": infer_condition(path),
        "sample": infer_sample(path),
        "cell": infer_cell(path),
        "file_name": path.name,
        "file_path": str(path),
    }


def discover_trackmate_spot_files(raw_data_dir: Path) -> pd.DataFrame:
    """Find TrackMate *_spots.csv files and return metadata table."""
    raw_data_dir = Path(raw_data_dir)
    spot_files = sorted(raw_data_dir.rglob("*_spots.csv"))

    metadata = pd.DataFrame(
        [parse_file_metadata(path) for path in spot_files]
    )

    if metadata.empty:
        return metadata

    return metadata.sort_values(
        ["protein", "sample", "cell", "condition", "file_name"],
        na_position="last",
    ).reset_index(drop=True)


def read_trackmate_spots(
    file_path: Path,
    frame_interval_s: float = 0.02,
) -> pd.DataFrame:
    """Read TrackMate spots CSV and return standardized columns.

    TrackMate CSV exports may contain service rows. This function standardizes
    column names and removes rows where required numeric fields cannot be parsed.
    """
    file_path = Path(file_path)

    raw = pd.read_csv(file_path, low_memory=False)
    raw.columns = [standardize_column_name(column) for column in raw.columns]

    rename_map = {
        column: COLUMN_ALIASES[column]
        for column in raw.columns
        if column in COLUMN_ALIASES
    }

    spots = raw.rename(columns=rename_map).copy()

    required_columns = ["track_id", "x", "y", "frame"]
    missing_columns = [
        column for column in required_columns if column not in spots.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {file_path.name}: "
            f"{missing_columns}. Available columns: {list(spots.columns)}"
        )

    numeric_columns = [
        "track_id",
        "x",
        "y",
        "frame",
        "t",
        "quality",
        "snr",
        "mean_intensity",
    ]

    for column in numeric_columns:
        if column in spots.columns:
            spots[column] = pd.to_numeric(spots[column], errors="coerce")

    spots = spots.dropna(subset=required_columns).copy()

    spots["track_id"] = spots["track_id"].astype(int)
    spots["frame"] = spots["frame"].astype(int)

    if "t" not in spots.columns or spots["t"].isna().all():
        spots["t"] = spots["frame"] * frame_interval_s

    return spots.sort_values(["track_id", "frame"]).reset_index(drop=True)


def save_table(dataframe: pd.DataFrame, output_path: Path) -> None:
    """Save DataFrame as CSV and create parent directory if needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(output_path, index=False)
