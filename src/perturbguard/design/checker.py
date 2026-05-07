from __future__ import annotations

import pandas as pd


def check_design(design: pd.DataFrame, min_controls_per_batch: int = 1, min_cells_per_perturbation: int = 50) -> pd.DataFrame:
    rows = []
    required = ["batch", "perturbation"]
    for column in required:
        if column not in design:
            rows.append(
                {
                    "check_id": f"required_{column}",
                    "status": "FAIL",
                    "message": f"Design is missing {column}.",
                    "suggested_fix": f"Add {column} to the design table.",
                }
            )
    if rows:
        return pd.DataFrame(rows)
    is_control = design["perturbation"].astype(str).str.lower().isin({"control", "ctrl", "ntc", "dmso"})
    for batch, group in design.groupby("batch", observed=True):
        n_controls = int(is_control.loc[group.index].sum())
        rows.append(
            {
                "check_id": "controls_per_batch",
                "status": "PASS" if n_controls >= min_controls_per_batch else "FAIL",
                "message": f"{batch} has {n_controls} planned control rows.",
                "suggested_fix": "Include controls in every batch.",
            }
        )
    if "n_cells_planned" in design:
        for perturbation, group in design.groupby("perturbation", observed=True):
            n_cells = int(group["n_cells_planned"].sum())
            rows.append(
                {
                    "check_id": "planned_cells_per_perturbation",
                    "status": "PASS" if n_cells >= min_cells_per_perturbation else "WARNING",
                    "message": f"{perturbation} has {n_cells} planned cells.",
                    "suggested_fix": "Increase planned cells for low-support perturbations.",
                }
            )
    return pd.DataFrame(rows)
