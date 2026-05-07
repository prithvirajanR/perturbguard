from __future__ import annotations

import pandas as pd
from anndata import AnnData


def _status_counts(sections: dict[str, pd.DataFrame]) -> str:
    counts: dict[str, int] = {}
    for df in sections.values():
        if "status" not in df:
            continue
        for status, count in df["status"].astype(str).value_counts().items():
            counts[status] = counts.get(status, 0) + int(count)
    if not counts:
        return "- No status-bearing checks were supplied."
    return "\n".join(f"- {status}: {count}" for status, count in sorted(counts.items()))


def build_dataset_card(
    adata: AnnData,
    sections: dict[str, pd.DataFrame],
    dataset_name: str = "dataset",
) -> str:
    perturbations = int(adata.obs["perturbation"].nunique()) if "perturbation" in adata.obs else 0
    controls = int(adata.obs["is_control"].sum()) if "is_control" in adata.obs else 0
    return f"""# Dataset Card: {dataset_name}

## Summary

- n_cells: {adata.n_obs}
- n_genes: {adata.n_vars}
- n_perturbations: {perturbations}
- n_controls: {controls}

## Audit Status Counts

{_status_counts(sections)}

## Supported Uses

- Dataset QC and benchmark split auditing.
- Perturbation generalization claim checking when a split is provided.
- Reporting known confounding, control, support, guide, and target-mapping risks.

## Known Limitations

- PASS/WARNING/FAIL thresholds are guardrails, not biological truth.
- Missing optional metadata reduces check coverage.
- Drug/pathway targets may require user-supplied mapping tables.
"""
