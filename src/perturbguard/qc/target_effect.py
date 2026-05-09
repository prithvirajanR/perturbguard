from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbguard.utils.matrix import mean_axis0
from perturbguard.qc.matched_controls import find_matched_controls
from perturbguard.qc.target_mapping import classify_target
from perturbguard.utils.controls import control_mask


def _expected_direction(group: pd.DataFrame) -> str:
    if "expected_direction" in group:
        values = group["expected_direction"].dropna().astype(str).str.lower()
        values = values[values.isin(["up", "down", "any"])]
        if not values.empty:
            return str(values.iloc[0])
    if "perturbation_type" not in group:
        return "down"
    perturbation_type = str(group["perturbation_type"].dropna().astype(str).iloc[0]).lower()
    if any(token in perturbation_type for token in ["crispra", "activation", "overexpression", "orf"]):
        return "up"
    if any(token in perturbation_type for token in ["crispri", "knockdown", "ko", "knockout", "sirna"]):
        return "down"
    return "any"


def _direction_status(logfc: float, expected_direction: str, fail_abs: float, warning_abs: float) -> tuple[str, str]:
    if expected_direction == "up":
        status = "PASS" if logfc >= warning_abs else "WARNING" if logfc > fail_abs else "FAIL"
        return status, "Target expression direction matches expected increase." if status == "PASS" else "Target expression does not show expected increase."
    if expected_direction == "down":
        status = "PASS" if logfc <= -warning_abs else "WARNING" if logfc < -fail_abs else "FAIL"
        return status, "Target expression direction matches expected decrease." if status == "PASS" else "Target expression does not show expected decrease."
    status = "PASS" if abs(logfc) >= warning_abs else "WARNING" if abs(logfc) > fail_abs else "FAIL"
    return status, "Target expression shows a detectable change." if status == "PASS" else "Target expression change is weak or absent."


def _target_genes(group: pd.DataFrame, fallback_target: str, gene_names: set[str]) -> tuple[list[str], str, str]:
    if "target_genes" in group:
        raw_values = group["target_genes"].dropna().astype(str)
        if not raw_values.empty:
            raw = raw_values.iloc[0]
            parsed = [
                value.strip()
                for chunk in raw.replace(";", ",").replace("|", ",").split(",")
                for value in [chunk]
                if value.strip()
            ]
            normalized = [gene if gene in gene_names else gene.upper() for gene in parsed]
            measured = [gene for gene in normalized if gene in gene_names]
            if measured:
                return measured, "multi_gene" if len(measured) > 1 else "gene", raw
    classification = classify_target(fallback_target, gene_names)
    if classification.target_type == "gene":
        return [classification.normalized_target], classification.target_type, fallback_target
    return [], classification.target_type, fallback_target


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
        target_genes, target_type, raw_target = _target_genes(group, target, gene_names)
        if not target_genes:
            target_classification = classify_target(target, gene_names)
            rows.append({"perturbation": pert, "target_gene": raw_target, "target_type": target_type, "n_target_genes": 0, "log2_fold_change": np.nan, "expected_direction": "down", "control_match_level": "", "status": "INSUFFICIENT_METADATA", "message": target_classification.message})
            continue
        idx = [var_index[gene] for gene in target_genes]
        matched_control_mask, match_level = find_matched_controls(obs, group.iloc[0])
        control_mean = mean_axis0(adata.X[matched_control_mask.values])
        group_mask = obs.index.isin(group.index)
        pert_mean = mean_axis0(adata.X[group_mask])[idx]
        logfc_values = np.log2((pert_mean + 1e-3) / (control_mean[idx] + 1e-3))
        logfc = float(np.nanmean(logfc_values))
        expected_direction = _expected_direction(group)
        status, message = _direction_status(
            logfc,
            expected_direction,
            logfc_fail_abs,
            logfc_warning_abs,
        )
        rows.append({"perturbation": pert, "target_gene": raw_target, "target_type": target_type, "n_target_genes": len(target_genes), "log2_fold_change": logfc, "expected_direction": expected_direction, "control_match_level": match_level, "status": status, "message": message})
    return pd.DataFrame(rows)
