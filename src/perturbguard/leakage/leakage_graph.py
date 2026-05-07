from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


def build_leakage_graph(adata: AnnData, split: pd.DataFrame) -> pd.DataFrame:
    obs = adata.obs.copy()
    split_series = align_split_to_obs(obs, split)
    valid_split, _ = validate_split_groups(split_series)
    if not valid_split:
        return pd.DataFrame(columns=["source", "target", "edge_type", "metadata_value"])
    rows = []
    for column in ["perturbation", "target_gene", "batch", "donor", "cell_type"]:
        if column not in obs:
            continue
        train_values = set(obs.loc[split_series.eq("train"), column].dropna().astype(str))
        test_values = set(obs.loc[split_series.isin(["val", "test"]), column].dropna().astype(str))
        for value in sorted(train_values & test_values):
            rows.append(
                {
                    "source": f"train:{value}",
                    "target": f"test:{value}",
                    "edge_type": f"{column}_overlap",
                    "metadata_value": value,
                }
            )
    return pd.DataFrame(rows, columns=["source", "target", "edge_type", "metadata_value"])
