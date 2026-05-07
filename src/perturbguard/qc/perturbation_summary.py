from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.utils.controls import control_mask


def summarize_perturbations(adata: AnnData, min_cells_warning: int = 50, min_cells_fail: int = 20) -> pd.DataFrame:
    rows = []
    obs = adata.obs
    controls = control_mask(obs)
    for pert, group in obs.groupby("perturbation", observed=True):
        n = len(group)
        status = "FAIL" if n < min_cells_fail else "WARNING" if n < min_cells_warning else "PASS"
        warnings = []
        if status != "PASS":
            warnings.append("low cell count")
        if "batch" in group and group["batch"].nunique() <= 1 and pert != "control":
            warnings.append("single-batch perturbation")
        matched_controls = int((obs.loc[controls, "batch"].isin(group.get("batch", pd.Series(dtype=str)).unique())).sum()) if "batch" in obs else int(controls.sum())
        rows.append({
            "perturbation": pert,
            "target_gene": ",".join(sorted(map(str, group.get("target_gene", pd.Series()).dropna().unique()))),
            "n_cells": n,
            "n_guides": group["guide_id"].nunique() if "guide_id" in group else 0,
            "n_batches": group["batch"].nunique() if "batch" in group else 0,
            "n_replicates": group["replicate"].nunique() if "replicate" in group else 0,
            "n_cell_types": group["cell_type"].nunique() if "cell_type" in group else 0,
            "n_controls_matched": matched_controls,
            "status": status,
            "warnings": "; ".join(warnings),
        })
    return pd.DataFrame(rows)
