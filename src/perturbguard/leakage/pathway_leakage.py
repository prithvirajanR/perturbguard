from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


def check_pathway_leakage(adata: AnnData, split: pd.DataFrame, pathways: dict[str, list[str]]) -> pd.DataFrame:
    obs = adata.obs.copy()
    split_series = align_split_to_obs(obs, split)
    valid_split, split_message = validate_split_groups(split_series)
    if not valid_split:
        return pd.DataFrame([{"pathway": "", "train_targets": "", "test_targets": "", "n_overlap": 0, "status": "INSUFFICIENT_METADATA", "message": split_message}])
    if "target_gene" not in obs:
        return pd.DataFrame(
            [
                {
                    "pathway": "",
                    "train_targets": "",
                    "test_targets": "",
                    "n_overlap": 0,
                    "status": "INSUFFICIENT_METADATA",
                }
            ]
        )
    train_targets = set(obs.loc[split_series.eq("train"), "target_gene"].dropna().astype(str))
    test_targets = set(obs.loc[split_series.isin(["val", "test"]), "target_gene"].dropna().astype(str))
    rows = []
    for pathway, genes in pathways.items():
        gene_set = set(map(str, genes))
        train = sorted(gene_set & train_targets)
        test = sorted(gene_set & test_targets)
        overlap = sorted(set(train) & set(test))
        rows.append(
            {
                "pathway": pathway,
                "train_targets": ",".join(train),
                "test_targets": ",".join(test),
                "n_overlap": len(overlap),
                "status": "WARNING" if train and test else "PASS",
            }
        )
    return pd.DataFrame(rows)
