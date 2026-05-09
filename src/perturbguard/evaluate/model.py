from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score

from perturbguard.utils.controls import control_mask
from perturbguard.utils.matrix import mean_axis0


def _status(value: float, pass_at: float = 0.8, warn_at: float = 0.5) -> str:
    if value >= pass_at:
        return "PASS"
    if value >= warn_at:
        return "WARNING"
    return "FAIL"


def _require_columns(table: pd.DataFrame, columns: set[str]) -> None:
    missing = sorted(columns - set(table.columns))
    if missing:
        raise ValueError(f"Predictions table is missing required columns: {missing}")


def _metrics(y_true: pd.Series, y_pred: pd.Series) -> dict[str, float]:
    if len(y_true) == 0:
        return {"accuracy": 0.0, "balanced_accuracy": 0.0, "macro_f1": 0.0}
    labels = sorted(set(y_true.astype(str)) | set(y_pred.astype(str)))
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)),
    }


def _prediction_genes(predictions: pd.DataFrame) -> list[str]:
    genes = []
    for column in predictions.columns:
        if column.startswith("pred_") and not column.startswith(("pred_low_", "pred_high_")):
            gene = column.removeprefix("pred_")
            if f"true_{gene}" in predictions:
                genes.append(gene)
    return genes


def _corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan
    value = spearmanr(x, y, nan_policy="omit").statistic
    return float(value) if value is not None else np.nan


def _effect_metrics(
    adata: AnnData,
    predictions: pd.DataFrame,
    genes: list[str],
    top_k_de: int,
) -> pd.DataFrame:
    obs = adata.obs
    if "perturbation" not in predictions:
        return pd.DataFrame(
            [
                {
                    "check_id": "perturbation_effect_recovery",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "Expression prediction rows need a perturbation column.",
                }
            ]
        )
    available_genes = [gene for gene in genes if gene in adata.var_names]
    if len(available_genes) < 2:
        return pd.DataFrame(
            [
                {
                    "check_id": "perturbation_effect_recovery",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "At least two predicted genes must match AnnData var names.",
                }
            ]
        )
    var_index = {str(gene): i for i, gene in enumerate(adata.var_names)}
    controls = control_mask(obs)
    if not controls.any():
        return pd.DataFrame(
            [
                {
                    "check_id": "perturbation_effect_recovery",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "No controls are available for perturbation-effect recovery.",
                }
            ]
        )
    control_mean = mean_axis0(adata.X[controls.values])
    rows = []
    for _, row in predictions.iterrows():
        perturbation = str(row["perturbation"])
        group_mask = obs["perturbation"].astype(str).eq(perturbation) & ~controls
        if group_mask.sum() == 0:
            rows.append(
                {
                    "check_id": "perturbation_effect_recovery",
                    "perturbation": perturbation,
                    "status": "INSUFFICIENT_METADATA",
                    "n_genes": len(available_genes),
                    "message": f"No observed perturbed cells found for {perturbation}.",
                }
            )
            continue
        idx = [var_index[gene] for gene in available_genes]
        observed_effect = mean_axis0(adata.X[group_mask.values])[idx] - control_mean[idx]
        predicted_effect = row[[f"pred_{gene}" for gene in available_genes]].astype(float).to_numpy() - control_mean[idx]
        if all(f"true_{gene}" in predictions for gene in available_genes):
            true_values = row[[f"true_{gene}" for gene in available_genes]].astype(float).to_numpy()
            observed_effect = true_values - control_mean[idx]
        k = min(top_k_de, len(available_genes))
        observed_top = set(np.argsort(np.abs(observed_effect))[-k:])
        predicted_top = set(np.argsort(np.abs(predicted_effect))[-k:])
        recovery = len(observed_top & predicted_top) / max(k, 1)
        corr = _corr(observed_effect, predicted_effect)
        direction = float(np.mean(np.sign(observed_effect) == np.sign(predicted_effect)))
        combined = np.nanmean([corr if not np.isnan(corr) else 0.0, recovery, direction])
        rows.append(
            {
                "check_id": "perturbation_effect_recovery",
                "perturbation": perturbation,
                "status": _status(float(combined), pass_at=0.7, warn_at=0.4),
                "n_genes": len(available_genes),
                "spearman_effect_corr": corr,
                "top_k_de_recovery": recovery,
                "direction_agreement": direction,
                "message": (
                    f"{perturbation}: effect rank correlation {corr:.2f}, "
                    f"top-{k} DE recovery {recovery:.2f}."
                ),
            }
        )
    return pd.DataFrame(rows)


def _pathway_metrics(
    adata: AnnData,
    predictions: pd.DataFrame,
    genes: list[str],
    gene_sets: dict[str, list[str]] | None,
) -> pd.DataFrame:
    if not gene_sets:
        return pd.DataFrame(
            [
                {
                    "check_id": "pathway_recovery",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "No gene sets were provided for pathway recovery.",
                }
            ]
        )
    controls = control_mask(adata.obs)
    if not controls.any() or "perturbation" not in predictions:
        return pd.DataFrame(
            [
                {
                    "check_id": "pathway_recovery",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "Controls and perturbation-level rows are required for pathway recovery.",
                }
            ]
        )
    control_mean = mean_axis0(adata.X[controls.values])
    var_index = {str(gene): i for i, gene in enumerate(adata.var_names)}
    rows = []
    for _, row in predictions.iterrows():
        perturbation = str(row["perturbation"])
        observed_scores = []
        predicted_scores = []
        for pathway, members in gene_sets.items():
            pathway_genes = [gene for gene in members if gene in genes and gene in var_index]
            if not pathway_genes:
                continue
            idx = [var_index[gene] for gene in pathway_genes]
            observed = row[[f"true_{gene}" for gene in pathway_genes]].astype(float).to_numpy()
            predicted = row[[f"pred_{gene}" for gene in pathway_genes]].astype(float).to_numpy()
            observed_scores.append(float(np.mean(observed - control_mean[idx])))
            predicted_scores.append(float(np.mean(predicted - control_mean[idx])))
        corr = _corr(np.asarray(observed_scores), np.asarray(predicted_scores))
        status = "INSUFFICIENT_METADATA" if np.isnan(corr) else _status(corr, pass_at=0.7, warn_at=0.4)
        rows.append(
            {
                "check_id": "pathway_recovery",
                "perturbation": perturbation,
                "status": status,
                "n_pathways": len(observed_scores),
                "pathway_score_correlation": corr,
                "message": f"{perturbation}: pathway score correlation {corr:.2f}.",
            }
        )
    return pd.DataFrame(rows)


def _interval_metrics(predictions: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    interval_genes = [
        gene for gene in genes if f"pred_low_{gene}" in predictions and f"pred_high_{gene}" in predictions
    ]
    if not interval_genes:
        return pd.DataFrame(
            [
                {
                    "check_id": "prediction_interval_calibration",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "No pred_low_/pred_high_ gene columns were provided.",
                }
            ]
        )
    lows = predictions[[f"pred_low_{gene}" for gene in interval_genes]].astype(float).to_numpy()
    highs = predictions[[f"pred_high_{gene}" for gene in interval_genes]].astype(float).to_numpy()
    truth = predictions[[f"true_{gene}" for gene in interval_genes]].astype(float).to_numpy()
    covered = (truth >= lows) & (truth <= highs)
    coverage = float(np.mean(covered))
    return pd.DataFrame(
        [
            {
                "check_id": "prediction_interval_calibration",
                "status": "PASS" if coverage >= 0.8 else "WARNING" if coverage >= 0.6 else "FAIL",
                "n_values": int(covered.size),
                "interval_coverage": coverage,
                "message": f"Prediction interval coverage is {coverage:.2f}.",
            }
        ]
    )


def evaluate_predictions(
    adata: AnnData,
    predictions: pd.DataFrame,
    group_columns: list[str] | None = None,
    min_group_size: int = 20,
    gene_sets: dict[str, list[str]] | None = None,
    top_k_de: int = 50,
) -> dict[str, pd.DataFrame]:
    """Evaluate perturbation predictions with audit-friendly status tables."""
    genes = _prediction_genes(predictions)
    if genes:
        return {
            "effect_recovery": _effect_metrics(adata, predictions, genes, top_k_de),
            "pathway_recovery": _pathway_metrics(adata, predictions, genes, gene_sets),
            "interval_calibration": _interval_metrics(predictions, genes),
        }

    _require_columns(predictions, {"cell_id", "y_true", "y_pred"})
    pred = predictions.copy()
    pred["cell_id"] = pred["cell_id"].astype(str)
    obs = adata.obs.copy()
    obs["_cell_id"] = obs.index.astype(str)
    merged = pred.merge(obs, left_on="cell_id", right_on="_cell_id", how="left", validate="one_to_one")
    missing_obs = int(merged["_cell_id"].isna().sum())
    metric_values = _metrics(merged["y_true"].astype(str), merged["y_pred"].astype(str))
    overall = pd.DataFrame(
        [
            {
                "check_id": "overall_model_performance",
                "status": _status(metric_values["balanced_accuracy"]),
                "n_predictions": int(len(merged)),
                "missing_cells": missing_obs,
                **metric_values,
                "message": (
                    f"Balanced accuracy {metric_values['balanced_accuracy']:.2f}; "
                    f"{missing_obs} predictions could not be matched to AnnData cells."
                ),
            }
        ]
    )

    groups = []
    for column in group_columns or ["batch", "cell_type", "perturbation"]:
        if column not in merged:
            continue
        for value, group in merged.groupby(column, observed=True, dropna=True):
            if len(group) < min_group_size:
                groups.append(
                    {
                        "check_id": "group_model_performance",
                        "group_column": column,
                        "group_value": str(value),
                        "status": "INSUFFICIENT_METADATA",
                        "n_predictions": int(len(group)),
                        "balanced_accuracy": 0.0,
                        "macro_f1": 0.0,
                        "message": f"{column}={value} has fewer than {min_group_size} predictions.",
                    }
                )
                continue
            values = _metrics(group["y_true"].astype(str), group["y_pred"].astype(str))
            groups.append(
                {
                    "check_id": "group_model_performance",
                    "group_column": column,
                    "group_value": str(value),
                    "status": _status(values["balanced_accuracy"]),
                    "n_predictions": int(len(group)),
                    **values,
                    "message": f"{column}={value} balanced accuracy {values['balanced_accuracy']:.2f}.",
                }
            )
    per_group = pd.DataFrame(groups)
    if per_group.empty:
        per_group = pd.DataFrame(
            [
                {
                    "check_id": "group_model_performance",
                    "status": "INSUFFICIENT_METADATA",
                    "message": "No requested group columns were available.",
                }
            ]
        )

    if "confidence" in merged:
        correct = merged["y_true"].astype(str).eq(merged["y_pred"].astype(str)).astype(float)
        confidence = pd.to_numeric(merged["confidence"], errors="coerce")
        valid = confidence.notna()
        if valid.sum() >= max(10, min_group_size):
            error = float(np.abs(confidence[valid] - correct[valid]).mean())
            calibration = pd.DataFrame(
                [
                    {
                        "check_id": "confidence_calibration",
                        "status": "PASS" if error <= 0.2 else "WARNING" if error <= 0.35 else "FAIL",
                        "n_predictions": int(valid.sum()),
                        "mean_absolute_calibration_error": error,
                        "message": f"Mean absolute confidence error {error:.2f}.",
                    }
                ]
            )
        else:
            calibration = pd.DataFrame(
                [
                    {
                        "check_id": "confidence_calibration",
                        "status": "INSUFFICIENT_METADATA",
                        "n_predictions": int(valid.sum()),
                        "message": "Not enough numeric confidence values to audit calibration.",
                    }
                ]
            )
    else:
        calibration = pd.DataFrame(
            [
                {
                    "check_id": "confidence_calibration",
                    "status": "WARNING",
                    "message": "Predictions table has no confidence column.",
                }
            ]
        )
    return {"overall": overall, "per_group": per_group, "calibration": calibration}
