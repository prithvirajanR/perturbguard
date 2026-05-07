from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.utils.controls import control_mask


def compare_datasets(datasets: dict[str, AnnData]) -> pd.DataFrame:
    rows = []
    for name, adata in datasets.items():
        rows.append(
            {
                "dataset": name,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_perturbations": int(adata.obs["perturbation"].nunique())
                if "perturbation" in adata.obs
                else 0,
                "n_batches": int(adata.obs["batch"].nunique()) if "batch" in adata.obs else 0,
                "n_controls": int(control_mask(adata.obs).sum()),
            }
        )
    return pd.DataFrame(rows)
