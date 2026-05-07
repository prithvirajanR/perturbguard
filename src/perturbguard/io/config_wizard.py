from __future__ import annotations

from collections.abc import Iterable

from anndata import AnnData


ALIASES = {
    "perturbation": ["perturbation", "condition", "gene", "drug", "treatment", "guide_target"],
    "is_control": ["is_control", "control", "is_ctrl"],
    "target_gene": ["target_gene", "target", "gene_target", "target_name"],
    "guide_id": ["guide_id", "guide", "grna", "sgRNA", "guide_name"],
    "batch": ["batch", "batch_id", "plate_batch"],
    "replicate": ["replicate", "rep", "bio_rep"],
    "cell_type": ["cell_type", "celltype", "cell_line"],
    "donor": ["donor", "donor_id"],
    "plate": ["plate", "plate_id", "well_plate"],
    "timepoint": ["timepoint", "time", "time_point"],
    "dose": ["dose", "concentration", "drug_dose"],
}


def _find_column(columns: Iterable[str], aliases: list[str]) -> str | None:
    by_lower = {column.lower(): column for column in columns}
    for alias in aliases:
        if alias.lower() in by_lower:
            return by_lower[alias.lower()]
    for column in columns:
        lowered = column.lower()
        if any(alias.lower() in lowered for alias in aliases):
            return column
    return None


def infer_config(adata: AnnData) -> dict:
    schema = {}
    for canonical, aliases in ALIASES.items():
        found = _find_column(adata.obs.columns.astype(str), aliases)
        if found and found != canonical:
            schema[canonical] = found
        elif found:
            schema[canonical] = canonical
    controls = {"labels": ["control", "ctrl", "ntc", "dmso", "vehicle"]}
    return {
        "schema": schema,
        "controls": controls,
        "notes": [
            "Review inferred mappings before using them for publication-grade audits.",
            "Columns are inferred from common perturbation dataset aliases.",
        ],
    }
