from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from anndata import AnnData


TARGET_CLASS_HINTS = {
    "hdac",
    "jak",
    "aurora kinase",
    "histone methyltransferase",
    "parp",
    "topoisomerase",
    "mek",
    "egfr",
    "vegfr",
    "cdk",
}


@dataclass(frozen=True)
class TargetClassification:
    raw_target: str
    target_type: str
    normalized_target: str
    message: str


def classify_target(target: object, gene_names: set[str]) -> TargetClassification:
    raw = "" if target is None else str(target).strip()
    if not raw or raw.lower() in {"nan", "none", "unknown"}:
        return TargetClassification(raw, "missing", "", "No target annotation is available.")
    if raw in gene_names or raw.upper() in gene_names:
        normalized = raw if raw in gene_names else raw.upper()
        return TargetClassification(raw, "gene", normalized, "Target maps directly to a measured gene.")
    lowered = raw.lower()
    if "/" in raw or "&" in raw or "pathway" in lowered or "signaling" in lowered:
        return TargetClassification(raw, "pathway_or_class", raw, "Target appears to describe a pathway or class, not a single gene.")
    if any(hint in lowered for hint in TARGET_CLASS_HINTS) or "," in raw:
        return TargetClassification(raw, "target_class", raw, "Target appears to describe a drug target class.")
    return TargetClassification(raw, "unmapped", raw, "Target does not map to a measured gene.")


def _mapping_lookup(target_map: pd.DataFrame | None) -> dict[str, dict[str, str]]:
    if target_map is None or target_map.empty or "perturbation" not in target_map or "target" not in target_map:
        return {}
    lookup = {}
    for _, row in target_map.iterrows():
        perturbation = str(row["perturbation"])
        lookup[perturbation] = {
            "target": str(row["target"]),
            "target_type": str(row.get("target_type", "")),
            "source": str(row.get("source", "user")),
        }
    return lookup


def audit_target_mapping(adata: AnnData, target_map: pd.DataFrame | None = None) -> pd.DataFrame:
    """Classify target annotations as gene, class/pathway, missing, or unmapped."""
    gene_names = set(adata.var_names.astype(str))
    lookup = _mapping_lookup(target_map)
    rows = []
    raw_targets: list[tuple[str, str, str]] = []
    if "target_gene" in adata.obs:
        raw_targets.extend((str(value), "target_gene", "") for value in adata.obs["target_gene"].dropna().unique())
    if "perturbation" in adata.obs:
        for perturbation in adata.obs["perturbation"].dropna().astype(str).unique():
            if perturbation in lookup:
                mapped = lookup[perturbation]
                raw_targets.append((mapped["target"], "target_map", mapped["source"]))
            elif perturbation not in {"control", "ctrl", "ntc", "DMSO", "dmso"}:
                raw_targets.append((perturbation, "perturbation", ""))
    seen = set()
    for raw, source_column, source in raw_targets:
        key = (raw, source_column, source)
        if key in seen:
            continue
        seen.add(key)
        classification = classify_target(raw, gene_names)
        status = "PASS" if classification.target_type == "gene" else "WARNING"
        if classification.target_type in {"missing", "unmapped"}:
            status = "INSUFFICIENT_METADATA"
        rows.append(
            {
                "check_id": "target_mapping",
                "status": status,
                "raw_target": classification.raw_target,
                "target_type": classification.target_type,
                "normalized_target": classification.normalized_target,
                "source_column": source_column,
                "source": source,
                "message": classification.message,
            }
        )
    if not rows:
        rows.append(
            {
                "check_id": "target_mapping",
                "status": "INSUFFICIENT_METADATA",
                "raw_target": "",
                "target_type": "missing",
                "normalized_target": "",
                "source_column": "",
                "source": "",
                "message": "No target annotations or perturbation labels were available.",
            }
        )
    return pd.DataFrame(rows)
