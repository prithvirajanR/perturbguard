from __future__ import annotations

import pandas as pd
from anndata import AnnData
from perturbguard.utils.controls import control_mask


def _status(n_cells: int, warning: int, fail: int) -> str:
    if n_cells < fail:
        return "FAIL"
    if n_cells < warning:
        return "WARNING"
    return "PASS"


def _row(entity: str, entity_type: str, n_cells: int, status: str, message: str) -> dict:
    return {
        "entity": entity,
        "entity_type": entity_type,
        "n_cells": int(n_cells),
        "status": status,
        "message": message,
    }


def check_cell_counts(
    adata: AnnData,
    min_cells_warning: int = 50,
    min_cells_fail: int = 20,
    min_cells_per_batch_warning: int = 10,
    min_cells_per_guide_warning: int = 10,
) -> pd.DataFrame:
    obs = adata.obs
    rows = []
    for perturbation, group in obs.groupby("perturbation", observed=True):
        n_cells = len(group)
        status = _status(n_cells, min_cells_warning, min_cells_fail)
        rows.append(
            _row(
                str(perturbation),
                "perturbation",
                n_cells,
                status,
                "Per-perturbation support is sufficient."
                if status == "PASS"
                else "Per-perturbation metrics may be unstable.",
            )
        )
    if "batch" in obs:
        for (perturbation, batch), group in obs.groupby(["perturbation", "batch"], observed=True):
            n_cells = len(group)
            status = "WARNING" if n_cells < min_cells_per_batch_warning else "PASS"
            rows.append(
                _row(
                    f"{perturbation}|{batch}",
                    "perturbation_batch",
                    n_cells,
                    status,
                    "Perturbation has enough cells in this batch."
                    if status == "PASS"
                    else "Perturbation has sparse batch-level support.",
                )
            )
    if "guide_id" in obs:
        for guide, group in obs.groupby("guide_id", observed=True):
            n_cells = len(group)
            status = "WARNING" if n_cells < min_cells_per_guide_warning else "PASS"
            rows.append(
                _row(
                    str(guide),
                    "guide",
                    n_cells,
                    status,
                    "Guide has enough cells."
                    if status == "PASS"
                    else "Guide-level effect estimates may be unstable.",
                )
            )
    if "is_control" in obs or "perturbation" in obs:
        n_controls = int(control_mask(obs).sum())
        status = _status(n_controls, min_cells_warning, min_cells_fail)
        rows.append(
            _row(
                "controls",
                "control",
                n_controls,
                status,
                "Control support is sufficient." if status == "PASS" else "Control support is sparse.",
            )
        )
    return pd.DataFrame(rows)
