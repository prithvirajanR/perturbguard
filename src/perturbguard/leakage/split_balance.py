from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


def _proportions(values: pd.Series) -> pd.Series:
    counts = values.astype(str).value_counts(normalize=True)
    return counts.sort_index()


def evaluate_split_balance(
    adata: AnnData,
    split: pd.DataFrame,
    variables: list[str] | None = None,
    warning: float = 0.25,
) -> pd.DataFrame:
    obs = adata.obs
    split_series = align_split_to_obs(obs, split)
    variables = variables or ["batch", "cell_type", "replicate"]
    valid_split, split_message = validate_split_groups(split_series)
    if not valid_split:
        return pd.DataFrame(
            [
                {
                    "variable": variable,
                    "max_proportion_difference": np.nan,
                    "status": "INSUFFICIENT_METADATA",
                    "message": split_message,
                }
                for variable in variables
            ]
        )
    rows = []
    for variable in variables:
        if variable not in obs:
            rows.append(
                {
                    "variable": variable,
                    "max_proportion_difference": np.nan,
                    "status": "INSUFFICIENT_METADATA",
                    "message": f"{variable} metadata is unavailable.",
                }
            )
            continue
        train = _proportions(obs.loc[split_series.eq("train"), variable])
        test = _proportions(obs.loc[split_series.isin(["val", "test"]), variable])
        levels = sorted(set(train.index) | set(test.index))
        max_diff = max(abs(float(train.get(level, 0.0)) - float(test.get(level, 0.0))) for level in levels)
        status = "WARNING" if max_diff > warning else "PASS"
        rows.append(
            {
                "variable": variable,
                "max_proportion_difference": max_diff,
                "status": status,
                "message": f"Maximum train/test proportion difference for {variable} is {max_diff:.2f}.",
            }
        )
    return pd.DataFrame(rows)
