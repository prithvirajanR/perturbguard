from __future__ import annotations

import pandas as pd
from anndata import AnnData

from perturbguard.io.schema import apply_schema


def _action(check_id: str, status: str, message: str) -> dict[str, str]:
    return {"check_id": check_id, "status": status, "message": message}


def repair_anndata(
    adata: AnnData,
    schema: dict[str, str] | None = None,
    control_labels: list[str] | None = None,
) -> tuple[AnnData, pd.DataFrame]:
    """Return a repaired in-memory AnnData plus a table of actions taken."""
    fixed = apply_schema(adata, schema or {}, control_labels).copy()
    actions: list[dict[str, str]] = []
    if not fixed.obs_names.is_unique:
        fixed.obs_names_make_unique()
        actions.append(_action("obs_names_unique", "WARNING", "Made duplicate cell IDs unique."))
    else:
        actions.append(_action("obs_names_unique", "PASS", "Cell IDs were already unique."))
    if not fixed.var_names.is_unique:
        fixed.var_names_make_unique()
        actions.append(_action("var_names_unique", "WARNING", "Made duplicate gene IDs unique."))
    else:
        actions.append(_action("var_names_unique", "PASS", "Gene IDs were already unique."))
    if "perturbation" in fixed.obs and "is_control" not in fixed.obs:
        labels = {label.lower() for label in (control_labels or ["control", "ctrl", "ntc", "dmso", "vehicle"])}
        fixed.obs["is_control"] = fixed.obs["perturbation"].astype(str).str.lower().isin(labels)
        actions.append(_action("is_control", "WARNING", "Inferred is_control from perturbation labels."))
    elif "is_control" in fixed.obs:
        fixed.obs["is_control"] = fixed.obs["is_control"].astype(bool)
        actions.append(_action("is_control", "PASS", "Canonical is_control column is present."))
    else:
        actions.append(_action("is_control", "FAIL", "Could not infer controls without perturbation labels."))
    for column in fixed.obs.columns:
        if fixed.obs[column].dtype == object:
            fixed.obs[column] = fixed.obs[column].astype("string")
    actions.append(_action("obs_dtypes", "PASS", "Coerced object observation columns to string dtype."))
    return fixed, pd.DataFrame(actions)
