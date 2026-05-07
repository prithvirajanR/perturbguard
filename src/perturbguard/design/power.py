from __future__ import annotations

import pandas as pd

from perturbguard.design.checker import check_design


def check_power(
    design: pd.DataFrame,
    min_cells_per_perturbation: int = 50,
    min_replicates: int = 2,
    min_batches: int = 2,
) -> pd.DataFrame:
    rows = check_design(
        design,
        min_controls_per_batch=1,
        min_cells_per_perturbation=min_cells_per_perturbation,
    ).to_dict("records")
    required = {"perturbation", "n_cells_planned"}
    if not required.issubset(design.columns):
        rows.append(
            {
                "check_id": "power_inputs",
                "status": "FAIL",
                "message": "Power checks require perturbation and n_cells_planned columns.",
                "suggested_fix": "Add planned cell counts per perturbation.",
            }
        )
        return pd.DataFrame(rows)
    for perturbation, group in design.groupby("perturbation", observed=True):
        n_cells = int(group["n_cells_planned"].sum())
        status = "PASS" if n_cells >= min_cells_per_perturbation else "WARNING"
        rows.append(
            {
                "check_id": "planned_cells_per_perturbation",
                "status": status,
                "message": f"{perturbation} has {n_cells} planned cells.",
                "suggested_fix": "Increase cells or downgrade claims for low-support perturbations.",
            }
        )
        if "replicate" in group:
            n_reps = int(group["replicate"].nunique())
            rows.append(
                {
                    "check_id": "replicate_support",
                    "status": "PASS" if n_reps >= min_replicates else "WARNING",
                    "message": f"{perturbation} spans {n_reps} replicate(s).",
                    "suggested_fix": "Use multiple biological/technical replicates per perturbation.",
                }
            )
        if "batch" in group:
            n_batches = int(group["batch"].nunique())
            rows.append(
                {
                    "check_id": "batch_support",
                    "status": "PASS" if n_batches >= min_batches else "WARNING",
                    "message": f"{perturbation} spans {n_batches} batch(es).",
                    "suggested_fix": "Spread perturbations across batches to reduce confounding.",
                }
            )
    return pd.DataFrame(rows)
