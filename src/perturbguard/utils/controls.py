from __future__ import annotations

import pandas as pd


TRUE_VALUES = {"true", "1", "yes", "y", "control", "ctrl", "ntc", "non-targeting", "dmso", "vehicle"}
FALSE_VALUES = {"false", "0", "no", "n", "nan", "none", ""}


def control_mask(obs: pd.DataFrame, column: str = "is_control") -> pd.Series:
    if column in obs:
        values = obs[column]
        if values.dtype == bool:
            return values.fillna(False)
        normalized = values.astype(str).str.strip().str.lower()
        return normalized.isin(TRUE_VALUES) & ~normalized.isin(FALSE_VALUES)
    if "perturbation" in obs:
        return obs["perturbation"].astype(str).str.strip().str.lower().isin(TRUE_VALUES)
    return pd.Series(False, index=obs.index)
