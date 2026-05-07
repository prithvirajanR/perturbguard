from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


def profile_h5ad(path: str | Path) -> pd.DataFrame:
    rows = []
    data_path = Path(path)
    try:
        adata = ad.read_h5ad(data_path, backed="r")
    except Exception as exc:
        return pd.DataFrame(
            [{"check_id": "backed_load", "status": "FAIL", "value": "", "message": str(exc)}]
        )
    try:
        rows.append({"check_id": "backed_load", "status": "PASS", "value": str(data_path), "message": "File opens in backed mode."})
        rows.append({"check_id": "n_cells", "status": "PASS", "value": int(adata.n_obs), "message": f"{adata.n_obs} cells."})
        rows.append({"check_id": "n_genes", "status": "PASS", "value": int(adata.n_vars), "message": f"{adata.n_vars} genes."})
        rows.append({"check_id": "x_storage", "status": "PASS", "value": type(adata.X).__name__, "message": "Backed expression storage type."})
        obs_mem = int(adata.obs.memory_usage(deep=True).sum())
        rows.append({"check_id": "obs_memory_bytes", "status": "PASS", "value": obs_mem, "message": "Observation metadata memory footprint."})
    finally:
        adata.file.close()
    return pd.DataFrame(rows)
