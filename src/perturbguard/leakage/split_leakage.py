from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


def _sets(obs: pd.DataFrame, split: pd.Series, column: str) -> tuple[set, set]:
    train = set(obs.loc[split.eq("train"), column].dropna().astype(str))
    test = set(obs.loc[split.isin(["val", "test"]), column].dropna().astype(str))
    return train, test


def check_split_leakage(adata: AnnData, split: pd.DataFrame, claim: str = "unseen_perturbation") -> pd.DataFrame:
    obs = adata.obs.copy()
    split_series = align_split_to_obs(obs, split)
    valid_split, split_message = validate_split_groups(split_series)
    if not valid_split:
        return pd.DataFrame(
            [
                {
                    "leakage_type": "split_validity",
                    "n_overlap": 0,
                    "overlap_values": "",
                    "status": "INSUFFICIENT_METADATA",
                    "message": split_message,
                }
            ]
        )
    rows = []
    checks = [("perturbation", "perturbation_overlap")]
    if claim in {"unseen_target_gene", "strict_unseen_combination_components"}:
        checks.append(("target_gene", "target_gene_overlap"))
    if claim == "unseen_cell_type":
        checks.append(("cell_type", "cell_type_overlap"))
    if claim == "unseen_donor":
        checks.append(("donor", "donor_overlap"))
    if claim == "unseen_batch":
        checks.append(("batch", "batch_overlap"))
    for col, leak_type in checks:
        if col not in obs:
            rows.append({"leakage_type": leak_type, "n_overlap": 0, "overlap_values": "", "status": "INSUFFICIENT_METADATA", "message": f"{col} metadata is unavailable."})
            continue
        train, test = _sets(obs, split_series, col)
        overlap = sorted(train & test)
        rows.append({"leakage_type": leak_type, "n_overlap": len(overlap), "overlap_values": ",".join(overlap[:20]), "status": "FAIL" if overlap else "PASS", "message": f"{len(overlap)} overlapping {col} values between train and test."})
    return pd.DataFrame(rows)
