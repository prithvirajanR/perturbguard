from __future__ import annotations

import pandas as pd
from anndata import AnnData

from perturbguard.leakage.split_leakage import check_split_leakage
from perturbguard.leakage.combination_leakage import check_combination_leakage
from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


CLAIM_TO_LEAKAGE = {
    "unseen_perturbation": "perturbation_overlap",
    "unseen_target_gene": "target_gene_overlap",
    "unseen_cell_type": "cell_type_overlap",
    "unseen_donor": "donor_overlap",
    "unseen_batch": "batch_overlap",
}


def check_claim(adata: AnnData, split: pd.DataFrame, claim: str) -> pd.DataFrame:
    valid_split, split_message = validate_split_groups(align_split_to_obs(adata.obs, split))
    if not valid_split:
        return pd.DataFrame(
            [{"claim": claim, "status": "INSUFFICIENT_METADATA", "message": split_message}]
        )
    leakage = check_split_leakage(adata, split, claim=claim)
    target = CLAIM_TO_LEAKAGE.get(claim)
    if claim in {"unseen_combinations", "strict_unseen_combination_components"}:
        combo = check_combination_leakage(adata, split)
        if combo["status"].eq("INSUFFICIENT_METADATA").all():
            status = "INSUFFICIENT_METADATA"
        elif claim == "strict_unseen_combination_components":
            if "n_components" not in combo or not combo["n_components"].ge(2).any():
                status = "INSUFFICIENT_METADATA"
            else:
                status = "SUPPORTED" if combo["status"].eq("PASS").all() else "UNSUPPORTED"
        else:
            status = "SUPPORTED" if not combo["status"].eq("FAIL").any() else "UNSUPPORTED"
        return pd.DataFrame([{"claim": claim, "status": status, "message": f"{claim} is {status.lower().replace('_', ' ')} under the provided split."}])
    if target is None:
        return pd.DataFrame([{"claim": claim, "status": "INSUFFICIENT_METADATA", "message": "Claim is not implemented in the MVP."}])
    row = leakage.loc[leakage["leakage_type"].eq(target)]
    if row.empty or row["status"].iloc[0] == "INSUFFICIENT_METADATA":
        status = "INSUFFICIENT_METADATA"
    elif row["status"].iloc[0] == "PASS":
        status = "SUPPORTED"
    else:
        status = "UNSUPPORTED"
    return pd.DataFrame([{"claim": claim, "status": status, "message": f"{claim} is {status.lower().replace('_', ' ')} under the provided split."}])
