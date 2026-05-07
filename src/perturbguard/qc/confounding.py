from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def _cramers_v(table: pd.DataFrame) -> tuple[float, float, np.ndarray]:
    chi2, p, _, expected = chi2_contingency(table)
    n = table.to_numpy().sum()
    denom = n * max(min(table.shape[0] - 1, table.shape[1] - 1), 1)
    return float(np.sqrt(chi2 / denom)) if denom else 0.0, float(p), expected


def detect_confounding(adata: AnnData, metadata_columns: list[str] | None = None, warning: float = 0.3, fail: float = 0.5) -> pd.DataFrame:
    obs = adata.obs
    metadata_columns = metadata_columns or ["batch", "replicate", "donor", "plate", "timepoint"]
    rows = []
    pvals = []
    for col in metadata_columns:
        if col not in obs:
            rows.append({"metadata_variable": col, "cramers_v": np.nan, "p_value": np.nan, "q_value": np.nan, "n_levels": 0, "status": "INSUFFICIENT_METADATA", "message": f"{col} metadata is unavailable."})
            continue
        table = pd.crosstab(obs["perturbation"], obs[col])
        if table.empty or table.shape[0] < 2 or table.shape[1] < 2:
            rows.append({"metadata_variable": col, "cramers_v": np.nan, "p_value": np.nan, "q_value": np.nan, "n_levels": table.shape[1], "status": "INSUFFICIENT_METADATA", "message": f"{col} has insufficient non-missing variation for confounding analysis."})
            continue
        v, p, expected = _cramers_v(table)
        pvals.append((len(rows), p))
        underpowered = (table.sum(axis=1).min() < 5) or (table.sum(axis=0).min() < 5) or (expected < 5).mean() > 0.2
        status = "INSUFFICIENT_METADATA" if underpowered else "FAIL" if v >= fail else "WARNING" if v >= warning else "PASS"
        rows.append({"metadata_variable": col, "cramers_v": v, "p_value": p, "q_value": p, "n_levels": table.shape[1], "status": status, "message": f"Perturbation association with {col}: Cramer's V={v:.2f}."})
    if pvals:
        qs = multipletests([p for _, p in pvals], method="fdr_bh")[1]
        for (idx, _), q in zip(pvals, qs, strict=False):
            rows[idx]["q_value"] = float(q)
    return pd.DataFrame(rows)


def per_perturbation_concentration(
    adata: AnnData,
    metadata_variable: str,
    warning: float = 0.7,
    fail: float = 0.9,
) -> pd.DataFrame:
    obs = adata.obs
    if metadata_variable not in obs:
        return pd.DataFrame(
            [
                {
                    "perturbation": "",
                    "metadata_variable": metadata_variable,
                    "dominant_group": "",
                    "dominant_fraction": np.nan,
                    "status": "INSUFFICIENT_METADATA",
                    "message": f"{metadata_variable} metadata is unavailable.",
                }
            ]
        )
    rows = []
    for perturbation, group in obs.groupby("perturbation", observed=True):
        counts = group[metadata_variable].value_counts(normalize=True)
        if counts.empty:
            rows.append(
                {
                    "perturbation": perturbation,
                    "metadata_variable": metadata_variable,
                    "dominant_group": "",
                    "dominant_fraction": np.nan,
                    "status": "INSUFFICIENT_METADATA",
                    "message": f"{metadata_variable} is missing for {perturbation}.",
                }
            )
            continue
        dominant_group = str(counts.index[0])
        dominant_fraction = float(counts.iloc[0])
        status = "FAIL" if dominant_fraction >= fail else "WARNING" if dominant_fraction >= warning else "PASS"
        rows.append(
            {
                "perturbation": perturbation,
                "metadata_variable": metadata_variable,
                "dominant_group": dominant_group,
                "dominant_fraction": dominant_fraction,
                "status": status,
                "message": f"{dominant_fraction:.0%} of {perturbation} cells are in {dominant_group}.",
            }
        )
    return pd.DataFrame(rows)
