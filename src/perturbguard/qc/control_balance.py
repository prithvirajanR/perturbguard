from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.utils.controls import control_mask


def check_control_balance(adata: AnnData, min_controls_per_group: int = 30) -> pd.DataFrame:
    obs = adata.obs
    if "is_control" not in obs and "perturbation" not in obs:
        return pd.DataFrame([{"grouping_variable": "is_control", "group_value": "", "n_controls": 0, "n_perturbed": len(obs), "control_fraction": 0.0, "status": "INSUFFICIENT_METADATA", "message": "Control labels are unavailable."}])
    control = control_mask(obs)
    rows = []
    for column in ["batch", "replicate", "cell_type", "timepoint"]:
        if column not in obs:
            rows.append({"grouping_variable": column, "group_value": "", "n_controls": 0, "n_perturbed": 0, "control_fraction": 0.0, "status": "INSUFFICIENT_METADATA", "message": f"{column} metadata is unavailable."})
            continue
        for value, idx in obs.groupby(column, observed=True).groups.items():
            mask = obs.index.isin(idx)
            n_controls = int((control & mask).sum())
            n_pert = int((~control & mask).sum())
            status = "PASS" if n_controls >= min_controls_per_group else "WARNING" if n_controls > 0 else "FAIL"
            rows.append({"grouping_variable": column, "group_value": value, "n_controls": n_controls, "n_perturbed": n_pert, "control_fraction": n_controls / max(n_controls + n_pert, 1), "status": status, "message": "Controls are available." if status == "PASS" else "Controls are sparse or absent for this group."})
    return pd.DataFrame(rows)
