from __future__ import annotations

import pandas as pd
from anndata import AnnData

from perturbguard.qc.metadata_shortcut import metadata_shortcut_baseline


def run_adversarial_checks(
    adata: AnnData,
    split: pd.DataFrame | None = None,
    min_cells_per_class: int = 10,
    n_splits: int = 5,
) -> pd.DataFrame:
    rows = []
    shortcut = metadata_shortcut_baseline(
        adata,
        min_cells_per_class=min_cells_per_class,
        n_splits=n_splits,
    ).iloc[0]
    rows.append(
        {
            "check_id": "metadata_predicts_perturbation",
            "status": shortcut["status"],
            "score": float(shortcut.get("balanced_accuracy", 0.0)),
            "message": shortcut["message"],
        }
    )
    if split is not None and "split" in split:
        split_counts = split["split"].astype(str).str.lower().value_counts(normalize=True)
        largest = float(split_counts.max()) if not split_counts.empty else 0.0
        rows.append(
            {
                "check_id": "split_label_dominance",
                "status": "WARNING" if largest > 0.9 else "PASS",
                "score": largest,
                "message": f"Largest split label accounts for {largest:.2%} of rows.",
            }
        )
    else:
        rows.append(
            {
                "check_id": "split_label_dominance",
                "status": "INSUFFICIENT_METADATA",
                "score": 0.0,
                "message": "No split table was provided.",
            }
        )
    return pd.DataFrame(rows)
