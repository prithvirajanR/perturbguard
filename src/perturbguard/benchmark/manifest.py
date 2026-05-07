from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

from perturbguard.claims.claim_checker import check_claim
from perturbguard.io.split_loader import load_split


def _row(check_id: str, status: str, message: str) -> dict[str, str]:
    return {"check_id": check_id, "status": status, "message": message}


def _resolve(base: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else base / path


def check_benchmark_manifest(path: str | Path) -> pd.DataFrame:
    manifest_path = Path(path)
    rows = []
    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except Exception as exc:
        return pd.DataFrame([_row("manifest_file", "FAIL", f"Could not read manifest: {exc}")])
    if not isinstance(manifest, dict):
        return pd.DataFrame([_row("manifest_file", "FAIL", "Benchmark manifest must be a mapping.")])
    rows.append(_row("manifest_file", "PASS", "Benchmark manifest loaded."))
    for field in ["dataset", "split", "claim", "model", "metrics"]:
        rows.append(
            _row(
                f"required_{field}",
                "PASS" if field in manifest else "FAIL",
                f"{field} is {'present' if field in manifest else 'missing'}.",
            )
        )
    if "dataset" not in manifest or "split" not in manifest or "claim" not in manifest:
        return pd.DataFrame(rows)
    data_path = _resolve(manifest_path.parent, manifest["dataset"])
    split_path = _resolve(manifest_path.parent, manifest["split"])
    rows.append(_row("dataset_path", "PASS" if data_path.exists() else "FAIL", str(data_path)))
    rows.append(_row("split_path", "PASS" if split_path.exists() else "FAIL", str(split_path)))
    if data_path.exists() and split_path.exists():
        try:
            claim = check_claim(ad.read_h5ad(data_path), load_split(split_path), str(manifest["claim"]))
            status = str(claim["status"].iloc[0])
            rows.append(_row("claim_support", status, str(claim["message"].iloc[0])))
        except Exception as exc:
            rows.append(_row("claim_support", "FAIL", f"Could not evaluate claim: {exc}"))
    return pd.DataFrame(rows)
