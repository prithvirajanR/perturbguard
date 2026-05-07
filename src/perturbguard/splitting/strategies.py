from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData


def _split_components(perturbation: str) -> tuple[str, ...]:
    normalized = perturbation.replace(";", "+").replace("|", "+").replace("-", "+")
    return tuple(part.strip() for part in normalized.split("+") if part.strip())


def _balanced_sample(
    obs: pd.DataFrame,
    columns: list[str],
    fraction: float,
    rng: np.random.Generator,
    available_index: pd.Index | None = None,
) -> pd.Index:
    source = obs.loc[available_index] if available_index is not None else obs
    chosen = []
    safe_columns = [column for column in columns if column in source]
    if not safe_columns:
        n = max(1, int(len(source) * fraction))
        return pd.Index(rng.choice(source.index.to_numpy(), size=min(n, len(source)), replace=False))
    strata = source[safe_columns].astype("string").fillna("<missing>").agg("||".join, axis=1)
    for _, stratum_index in strata.groupby(strata, observed=True).groups.items():
        stratum_index = pd.Index(stratum_index)
        if len(stratum_index) < 2:
            continue
        n = int(round(len(stratum_index) * fraction))
        n = min(max(1, n), len(stratum_index) - 1)
        chosen.extend(rng.choice(stratum_index.to_numpy(), size=n, replace=False).tolist())
    if not chosen:
        n = max(1, int(len(source) * fraction))
        chosen = rng.choice(source.index.to_numpy(), size=min(n, len(source)), replace=False).tolist()
    return pd.Index(chosen)


def generate_split(
    adata: AnnData,
    strategy: str = "random",
    test_size: float = 0.2,
    val_size: float = 0.0,
    random_state: int = 0,
    metadata_column: str | None = None,
    balance_columns: list[str] | None = None,
) -> pd.DataFrame:
    strategy = strategy.replace("-", "_")
    rng = np.random.default_rng(random_state)
    obs = adata.obs.copy()
    if not obs.index.is_unique:
        raise ValueError("Split generation requires unique observation/cell IDs.")
    split = pd.Series("train", index=obs.index, name="split")
    if strategy == "random":
        order = rng.permutation(len(obs))
        n_test = max(1, int(len(obs) * test_size))
        split.iloc[order[:n_test]] = "test"
    elif strategy == "balanced_random":
        columns = ["perturbation"] + list(balance_columns or [])
        test_idx = _balanced_sample(obs, columns, test_size, rng)
        split.loc[test_idx] = "test"
    elif strategy == "leave_target_gene_out":
        targets = sorted([t for t in obs["target_gene"].dropna().astype(str).unique() if t.lower() != "control"])
        rng.shuffle(targets)
        n_targets = max(1, int(np.ceil(len(targets) * test_size)))
        holdout = set(targets[:n_targets])
        split.loc[obs["target_gene"].astype(str).isin(holdout)] = "test"
    elif strategy == "leave_perturbation_out":
        perturbations = sorted([p for p in obs["perturbation"].dropna().astype(str).unique() if p.lower() != "control"])
        rng.shuffle(perturbations)
        holdout = set(perturbations[: max(1, int(np.ceil(len(perturbations) * test_size)))])
        split.loc[obs["perturbation"].astype(str).isin(holdout)] = "test"
    elif strategy == "leave_metadata_out":
        if not metadata_column or metadata_column not in obs:
            raise ValueError("leave_metadata_out requires an existing metadata_column.")
        values = sorted(obs[metadata_column].dropna().astype(str).unique())
        rng.shuffle(values)
        holdout = set(values[: max(1, int(np.ceil(len(values) * test_size)))])
        split.loc[obs[metadata_column].astype(str).isin(holdout)] = "test"
    elif strategy == "strict_unseen_combination":
        if "target_gene" not in obs:
            raise ValueError("strict_unseen_combination requires target_gene metadata.")
        targets = sorted([t for t in obs["target_gene"].dropna().astype(str).unique() if t.lower() != "control"])
        rng.shuffle(targets)
        holdout = set(targets[: max(1, int(np.ceil(len(targets) * test_size)))])
        split.loc[obs["target_gene"].astype(str).isin(holdout)] = "test"
    elif strategy == "seen_component_unseen_combination":
        perturbations = sorted(obs["perturbation"].dropna().astype(str).unique())
        combo_components = {
            p: _split_components(p)
            for p in perturbations
            if len(_split_components(p)) >= 2
        }
        combos = sorted(combo_components)
        if not combos:
            raise ValueError("seen_component_unseen_combination requires at least one combination perturbation.")
        rng.shuffle(combos)
        n_holdout = max(1, int(np.ceil(len(combos) * test_size)))
        holdout: set[str] = set()
        available = set(perturbations)
        for combo in combos:
            candidate_holdout = holdout | {combo}
            train_perturbations = available - candidate_holdout
            train_components = {component for pert in train_perturbations for component in _split_components(pert)}
            if set(combo_components[combo]).issubset(train_components):
                holdout.add(combo)
            if len(holdout) >= n_holdout:
                break
        if not holdout:
            raise ValueError("No seen-component unseen-combination split is possible for this dataset.")
        split.loc[obs["perturbation"].astype(str).isin(holdout)] = "test"
    else:
        raise ValueError(f"Unsupported split strategy: {strategy}")
    if val_size:
        train_idx = split[split.eq("train")].index.to_numpy()
        n_val = int(len(obs) * val_size)
        val_idx = rng.choice(train_idx, size=min(n_val, len(train_idx)), replace=False)
        split.loc[val_idx] = "val"
    return pd.DataFrame({"cell_id": obs.index.astype(str), "split": split.values})
