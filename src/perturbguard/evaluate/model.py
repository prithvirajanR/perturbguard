from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score


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


def evaluate_predictions(
    adata: AnnData,
    predictions: pd.DataFrame,
    group_columns: list[str] | None = None,
    min_group_size: int = 20,
) -> dict[str, pd.DataFrame]:
    """Evaluate perturbation predictions with audit-friendly status tables."""
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
