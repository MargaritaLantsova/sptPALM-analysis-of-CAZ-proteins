"""Hidden Markov Model analysis utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd

try:
    from hmmlearn.hmm import GaussianHMM
except ImportError:  # pragma: no cover
    GaussianHMM = None


STATE_LABELS = {
    0: "slow",
    1: "intermediate",
    2: "fast",
}


def prepare_hmm_sequences(
    links: pd.DataFrame,
    min_links_per_track: int = 3,
    observation_column: str = "log10_velocity",
) -> tuple[np.ndarray | None, list[int], list[int]]:
    """Prepare concatenated HMM observation matrix and sequence lengths."""
    observations = []
    lengths = []
    index_order = []

    group_columns = [
        "condition",
        "sample",
        "cell",
        "file_name",
        "track_id",
    ]

    for _, group in links.groupby(group_columns, sort=False):
        group = group.sort_values("start_frame")

        if len(group) < min_links_per_track:
            continue

        values = group[observation_column].to_numpy(dtype=float)

        if not np.all(np.isfinite(values)):
            continue

        observations.append(values.reshape(-1, 1))
        lengths.append(len(values))
        index_order.extend(group.index.tolist())

    if not observations:
        return None, [], []

    return np.vstack(observations), lengths, index_order


def order_hmm_states_by_mean_velocity(
    result_table: pd.DataFrame,
    state_column: str = "hmm_state_raw",
    value_column: str = "log10_velocity",
) -> tuple[pd.DataFrame, dict[int, int]]:
    """Reorder HMM states from slow to fast based on mean velocity."""
    result_table = result_table.copy()

    state_order = (
        result_table.groupby(state_column)[value_column]
        .mean()
        .sort_values()
        .index
    )

    state_map = {
        old_state: new_state
        for new_state, old_state in enumerate(state_order)
    }

    result_table["hmm_state"] = result_table[state_column].map(state_map)
    result_table["hmm_state_label"] = result_table["hmm_state"].map(
        STATE_LABELS
    )

    return result_table, state_map


def fit_hmm_for_protein(
    links: pd.DataFrame,
    protein: str,
    n_states: int = 3,
    min_links_per_track: int = 3,
    random_state: int = 42,
    n_iter: int = 300,
) -> tuple[GaussianHMM, pd.DataFrame, dict[int, int]]:
    """Fit a Gaussian HMM to one protein and return classified links."""
    if GaussianHMM is None:
        raise ImportError("hmmlearn is required. Install with: pip install hmmlearn")

    protein_data = links[links["protein"] == protein].copy()
    protein_data["log10_velocity"] = np.log10(protein_data["velocity_um_s"])

    observations, lengths, index_order = prepare_hmm_sequences(
        protein_data,
        min_links_per_track=min_links_per_track,
    )

    if observations is None:
        return None, pd.DataFrame(), {}

    model = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        n_iter=n_iter,
        random_state=random_state,
    )

    model.fit(observations, lengths)
    states_raw = model.predict(observations, lengths)

    result = protein_data.loc[index_order].copy()
    result["hmm_state_raw"] = states_raw

    result, state_map = order_hmm_states_by_mean_velocity(result)

    return model, result, state_map


def reorder_transition_matrix(
    model: GaussianHMM,
    state_map: dict[int, int],
) -> np.ndarray:
    """Reorder transition matrix according to slow-to-fast state map."""
    inverse_map = {
        new_state: old_state
        for old_state, new_state in state_map.items()
    }

    ordered_old_states = [
        inverse_map[index]
        for index in range(len(inverse_map))
    ]

    transition_matrix = model.transmat_
    return transition_matrix[np.ix_(ordered_old_states, ordered_old_states)]


def transition_matrix_to_table(
    transition_matrix: np.ndarray,
    protein: str,
) -> pd.DataFrame:
    """Convert HMM transition matrix to long-format table."""
    rows = []

    for from_state in range(transition_matrix.shape[0]):
        for to_state in range(transition_matrix.shape[1]):
            rows.append(
                {
                    "protein": protein,
                    "from_state": from_state,
                    "from_state_label": STATE_LABELS.get(
                        from_state,
                        f"state_{from_state}",
                    ),
                    "to_state": to_state,
                    "to_state_label": STATE_LABELS.get(
                        to_state,
                        f"state_{to_state}",
                    ),
                    "transition_probability": transition_matrix[
                        from_state,
                        to_state,
                    ],
                }
            )

    return pd.DataFrame(rows)


def summarize_hmm_states_by_cell(hmm_links: pd.DataFrame) -> pd.DataFrame:
    """Create per-cell HMM state summary."""
    rows = []

    for keys, group in hmm_links.groupby(
        ["protein", "condition", "sample", "cell", "file_name"],
        observed=True,
    ):
        protein, condition, sample, cell, file_name = keys
        n_links = len(group)
        state_counts = group["hmm_state_label"].value_counts()

        rows.append(
            {
                "protein": protein,
                "condition": condition,
                "sample": sample,
                "cell": cell,
                "file_name": file_name,
                "n_links": n_links,
                "slow_fraction": state_counts.get("slow", 0) / n_links,
                "intermediate_fraction": (
                    state_counts.get("intermediate", 0) / n_links
                ),
                "fast_fraction": state_counts.get("fast", 0) / n_links,
                "slow_percent": 100 * state_counts.get("slow", 0) / n_links,
                "intermediate_percent": (
                    100 * state_counts.get("intermediate", 0) / n_links
                ),
                "fast_percent": 100 * state_counts.get("fast", 0) / n_links,
            }
        )

    return pd.DataFrame(rows)
