from __future__ import annotations

import pandas as pd
from perturbguard.utils.controls import control_mask


def _control_mask(obs: pd.DataFrame) -> pd.Series:
    return control_mask(obs)


def find_matched_controls(obs: pd.DataFrame, perturbed_row: pd.Series) -> tuple[pd.Series, str]:
    controls = _control_mask(obs)
    candidates = controls.copy()
    levels = [
        ("same_batch_cell_type_timepoint", ["batch", "cell_type", "timepoint"]),
        ("same_batch_cell_type", ["batch", "cell_type"]),
        ("same_batch", ["batch"]),
        ("same_cell_type_timepoint", ["cell_type", "timepoint"]),
        ("same_cell_type", ["cell_type"]),
    ]
    for level, columns in levels:
        if not all(column in obs and column in perturbed_row.index for column in columns):
            continue
        mask = candidates.copy()
        for column in columns:
            mask &= obs[column].eq(perturbed_row[column])
        if mask.any():
            return mask, level
    return controls, "global_controls"
