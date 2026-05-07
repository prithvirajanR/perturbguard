from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbguard.utils.matrix import mean_axis0
from perturbguard.utils.controls import control_mask


def check_guide_consistency(adata: AnnData, corr_warning: float = 0.2, corr_fail: float = 0.0) -> pd.DataFrame:
    obs = adata.obs
    if not {"guide_id", "target_gene", "is_control"}.issubset(obs.columns):
        return pd.DataFrame([{"target_gene": "", "n_guides": 0, "median_correlation": np.nan, "min_correlation": np.nan, "status": "INSUFFICIENT_METADATA", "message": "Guide metadata is unavailable."}])
    controls = control_mask(obs).values
    if not controls.any():
        return pd.DataFrame([{"target_gene": "", "n_guides": 0, "n_pairs": 0, "median_correlation": np.nan, "min_correlation": np.nan, "status": "INSUFFICIENT_METADATA", "message": "No control cells are available for guide consistency checks."}])
    baseline = mean_axis0(adata.X[controls])
    rows = []
    for gene, gene_obs in obs.loc[~control_mask(obs)].groupby("target_gene", observed=True):
        guides = list(gene_obs["guide_id"].unique())
        if len(guides) < 2:
            continue
        effects = []
        for guide in guides:
            mask = obs["guide_id"].eq(guide).values
            effects.append(mean_axis0(adata.X[mask]) - baseline)
        cors = []
        for i in range(len(effects)):
            for j in range(i + 1, len(effects)):
                if np.std(effects[i]) == 0 or np.std(effects[j]) == 0:
                    cors.append(np.nan)
                    continue
                cors.append(float(np.corrcoef(effects[i], effects[j])[0, 1]))
        finite_cors = [corr for corr in cors if np.isfinite(corr)]
        if not finite_cors:
            rows.append({"target_gene": gene, "n_guides": len(guides), "n_pairs": len(cors), "median_correlation": np.nan, "min_correlation": np.nan, "status": "INSUFFICIENT_METADATA", "message": "Guide effect correlations are undefined, likely because effect vectors are constant."})
            continue
        med = float(np.median(finite_cors))
        mn = float(np.min(finite_cors))
        status = "FAIL" if med <= corr_fail else "WARNING" if med <= corr_warning else "PASS"
        rows.append({"target_gene": gene, "n_guides": len(guides), "n_pairs": len(cors), "median_correlation": med, "min_correlation": mn, "status": status, "message": "Guide effects are consistent." if status == "PASS" else "Guide effects are weakly correlated."})
    return pd.DataFrame(rows)
