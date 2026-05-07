from __future__ import annotations

from pathlib import Path

import pandas as pd


VALID_SPLITS = {"train", "val", "validation", "test"}


def normalize_split_labels(split: pd.DataFrame) -> pd.DataFrame:
    normalized = split.copy()
    if "split" not in normalized:
        raise ValueError("Split table must contain a 'split' column.")
    normalized["split"] = normalized["split"].astype(str).str.lower().replace({"validation": "val"})
    invalid = sorted(set(normalized["split"]) - {"train", "val", "test"})
    if invalid:
        raise ValueError(f"Invalid split labels: {invalid}")
    return normalized


def load_split(path: str | Path) -> pd.DataFrame:
    split = pd.read_csv(path)
    if "split" not in split:
        raise ValueError("Split CSV must contain a 'split' column.")
    if "cell_id" not in split and "perturbation" not in split:
        raise ValueError("Split CSV must contain either 'cell_id' or 'perturbation'.")
    split = normalize_split_labels(split)
    if "cell_id" in split:
        split["cell_id"] = split["cell_id"].astype(str)
    return split


def align_split_to_obs(obs: pd.DataFrame, split: pd.DataFrame) -> pd.Series:
    split = normalize_split_labels(split)
    if "cell_id" in split:
        by_cell = (
            split.assign(cell_id=split["cell_id"].astype(str))
            .drop_duplicates("cell_id", keep="last")
            .set_index("cell_id")["split"]
        )
        aligned = pd.Series(obs.index.astype(str), index=obs.index).map(by_cell)
    elif "perturbation" in split:
        if "perturbation" not in obs:
            raise ValueError("Perturbation-level split requires canonical 'perturbation' metadata.")
        by_perturbation = (
            split.assign(perturbation=split["perturbation"].astype(str))
            .drop_duplicates("perturbation", keep="last")
            .set_index("perturbation")["split"]
        )
        aligned = obs["perturbation"].astype(str).map(by_perturbation)
    else:
        raise ValueError("Split table must contain either cell_id or perturbation.")
    return aligned.fillna("train")


def validate_split_groups(split_labels: pd.Series) -> tuple[bool, str]:
    normalized = normalize_split_labels(pd.DataFrame({"split": split_labels}))["split"]
    has_train = normalized.eq("train").any()
    has_eval = normalized.isin(["val", "test"]).any()
    if not has_train and not has_eval:
        return False, "Split contains neither train nor evaluation cells."
    if not has_train:
        return False, "Split contains no train cells."
    if not has_eval:
        return False, "Split contains no evaluation cells."
    return True, "Split contains both train and evaluation cells."
