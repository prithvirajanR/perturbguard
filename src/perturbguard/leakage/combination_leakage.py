from __future__ import annotations

import re

import pandas as pd
from anndata import AnnData
from perturbguard.io.split_loader import align_split_to_obs, validate_split_groups


def parse_components(label: str, separators: list[str] | None = None) -> set[str]:
    separators = separators or ["+", "|", ";", "-"]
    pattern = "|".join(re.escape(s) for s in separators)
    return {p.strip() for p in re.split(pattern, str(label)) if p.strip() and p.strip().lower() not in {"control", "ctrl", "ntc", "dmso"}}


def categorize_combination_leakage(train_perturbations: list[str], test_perturbations: list[str]) -> pd.DataFrame:
    train_set = set(map(str, train_perturbations))
    train_components = set().union(*(parse_components(p) for p in train_set)) if train_set else set()
    rows = []
    for pert in map(str, test_perturbations):
        comps = parse_components(pert)
        if pert in train_set:
            category = "EXACT_LEAKAGE"
        elif comps and comps <= train_components:
            category = "RECOMBINATION_OF_SEEN_COMPONENTS"
        elif comps and comps & train_components:
            category = "PARTIAL_COMPONENT_HOLDOUT"
        elif comps:
            category = "STRICT_UNSEEN_COMBINATION"
        else:
            category = "UNKNOWN"
        rows.append({"perturbation": pert, "components": "+".join(sorted(comps)), "n_components": len(comps), "category": category})
    return pd.DataFrame(rows)


def check_combination_leakage(adata: AnnData, split: pd.DataFrame) -> pd.DataFrame:
    obs = adata.obs.copy()
    split_series = align_split_to_obs(obs, split)
    valid_split, split_message = validate_split_groups(split_series)
    if not valid_split:
        return pd.DataFrame(
            [
                {
                    "perturbation": "",
                    "components": "",
                    "category": "UNKNOWN",
                    "status": "INSUFFICIENT_METADATA",
                    "message": split_message,
                }
            ]
        )
    if "perturbation" not in obs:
        return pd.DataFrame(
            [
                {
                    "perturbation": "",
                    "components": "",
                    "category": "UNKNOWN",
                    "status": "INSUFFICIENT_METADATA",
                }
            ]
        )
    train = sorted(obs.loc[split_series.eq("train"), "perturbation"].astype(str).unique())
    test = sorted(obs.loc[split_series.isin(["val", "test"]), "perturbation"].astype(str).unique())
    results = categorize_combination_leakage(train, test)
    status_map = {
        "EXACT_LEAKAGE": "FAIL",
        "RECOMBINATION_OF_SEEN_COMPONENTS": "WARNING",
        "PARTIAL_COMPONENT_HOLDOUT": "WARNING",
        "STRICT_UNSEEN_COMBINATION": "PASS",
        "UNKNOWN": "INSUFFICIENT_METADATA",
    }
    results["status"] = results["category"].map(status_map)
    return results
