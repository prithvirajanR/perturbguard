from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbguard.utils.matrix import mean_axis0
from perturbguard.qc.matched_controls import find_matched_controls
from perturbguard.qc.target_mapping import classify_target
from perturbguard.utils.controls import control_mask


def check_target_effect(adata: AnnData, logfc_fail_abs: float = 0.05, logfc_warning_abs: float = 0.25) -> pd.DataFrame:
    obs = adata.obs
    try:
        if adata.X.min() < 0:
            return pd.DataFrame(
                [
                    {
                        "perturbation": "",
                        "target_gene": "",
                        "target_type": "unknown",
                        "log2_fold_change": np.nan,
                        "expected_direction": "down",
                        "control_match_level": "none",
                        "status": "INSUFFICIENT_METADATA",
                        "message": "Target-effect log fold change requires nonnegative expression values.",
                    }
                ]
            )
    except ValueError:
        pass
    controls = control_mask(obs)
    if not controls.any():
        return pd.DataFrame(
            [
                {
                    "perturbation": "",
                    "target_gene": "",
                    "log2_fold_change": np.nan,
                    "expected_direction": "down",
                    "control_match_level": "none",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "No control cells are available for target-effect checks.",
                }
            ]
        )
    rows = []
    var_index = {str(g): i for i, g in enumerate(adata.var_names)}
    gene_names = set(map(str, adata.var_names))
    for pert, group in obs.loc[~controls].groupby("perturbation", observed=True):
        target = str(group["target_gene"].iloc[0]) if "target_gene" in group else pert
        target_classification = classify_target(target, gene_names)
        if target_classification.target_type != "gene":
            rows.append({"perturbation": pert, "target_gene": target, "target_type": target_classification.target_type, "log2_fold_change": np.nan, "expected_direction": "down", "control_match_level": "", "status": "INSUFFICIENT_METADATA", "message": target_classification.message})
            continue
        idx = var_index[target_classification.normalized_target]
        matched_control_mask, match_level = find_matched_controls(obs, group.iloc[0])
        control_mean = mean_axis0(adata.X[matched_control_mask.values])
        group_mask = obs.index.isin(group.index)
        pert_mean = mean_axis0(adata.X[group_mask])[idx]
        logfc = float(np.log2((pert_mean + 1e-3) / (control_mean[idx] + 1e-3)))
        status = "PASS" if logfc <= -logfc_warning_abs else "WARNING" if logfc < -logfc_fail_abs else "FAIL"
        rows.append({"perturbation": pert, "target_gene": target, "target_type": target_classification.target_type, "log2_fold_change": logfc, "expected_direction": "down", "control_match_level": match_level, "status": status, "message": "Target expression direction matches expectation." if status == "PASS" else "Target expression does not show expected decrease."})
    return pd.DataFrame(rows)
