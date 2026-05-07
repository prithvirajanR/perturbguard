from __future__ import annotations

import pandas as pd
from anndata import AnnData
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


def metadata_shortcut_baseline(adata: AnnData, min_cells_per_class: int = 10, n_splits: int = 5, random_state: int = 0) -> pd.DataFrame:
    obs = adata.obs.copy()
    features = [c for c in ["batch", "replicate", "donor", "plate", "timepoint"] if c in obs]
    if not features:
        return pd.DataFrame([{"features_used": "", "n_classes": 0, "n_cells": len(obs), "balanced_accuracy": 0.0, "macro_f1": 0.0, "chance_level": 0.0, "status": "INSUFFICIENT_METADATA", "message": "No technical metadata columns available."}])
    counts = obs["perturbation"].value_counts()
    keep = obs["perturbation"].isin(counts[counts >= min_cells_per_class].index)
    obs = obs.loc[keep]
    min_class = int(obs["perturbation"].value_counts().min()) if len(obs) else 0
    folds = min(n_splits, min_class)
    if folds < 2 or obs["perturbation"].nunique() < 2:
        status, bal, f1 = "INSUFFICIENT_METADATA", 0.0, 0.0
    else:
        enc_kwargs = {"handle_unknown": "ignore"}
        try:
            encoder = OneHotEncoder(sparse_output=True, **enc_kwargs)
        except TypeError:
            encoder = OneHotEncoder(sparse=True, **enc_kwargs)
        model = Pipeline([
            ("prep", ColumnTransformer([("cat", encoder, features)])),
            ("clf", LogisticRegression(max_iter=1000, class_weight="balanced")),
        ])
        cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random_state)
        pred = cross_val_predict(model, obs[features], obs["perturbation"], cv=cv)
        bal = float(balanced_accuracy_score(obs["perturbation"], pred))
        f1 = float(f1_score(obs["perturbation"], pred, average="macro"))
        status = "FAIL" if bal >= 0.6 else "WARNING" if bal >= 0.4 else "PASS"
    chance = 1 / max(obs["perturbation"].nunique(), 1)
    return pd.DataFrame([{"features_used": ",".join(features), "n_classes": obs["perturbation"].nunique(), "n_cells": len(obs), "balanced_accuracy": bal, "macro_f1": f1, "chance_level": chance, "status": status, "message": f"Metadata-only balanced accuracy {bal:.2f} vs chance {chance:.2f}."}])
