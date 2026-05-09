from __future__ import annotations

import pandas as pd


STATUS_RANK = {
    "FAIL": 4,
    "UNSUPPORTED": 4,
    "WARNING": 3,
    "INSUFFICIENT_METADATA": 2,
    "PASS": 1,
    "SUPPORTED": 1,
}


def _worst_status(table: pd.DataFrame | None) -> str:
    if table is None or "status" not in table or table.empty:
        return "INSUFFICIENT_METADATA"
    statuses = table["status"].dropna().astype(str)
    if statuses.empty:
        return "INSUFFICIENT_METADATA"
    return max(statuses, key=lambda status: STATUS_RANK.get(status, 0))


def _message(table: pd.DataFrame | None, fallback: str) -> str:
    if table is None or table.empty or "message" not in table:
        return fallback
    flagged = table.loc[table["status"].astype(str).isin(["FAIL", "WARNING", "UNSUPPORTED"])]
    source = flagged if not flagged.empty else table
    return str(source["message"].iloc[0])


def build_preflight_checks(sections: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = [
        {
            "check_id": "metadata_confounding",
            "question": "Are batch, replicate, donor, plate, or timepoint confounded with perturbation labels?",
            "source_section": "confounding, metadata_concentration",
            "status": max(
                _worst_status(sections.get("confounding")),
                _worst_status(sections.get("metadata_concentration")),
                key=lambda status: STATUS_RANK.get(status, 0),
            ),
            "message": _message(
                sections.get("confounding"),
                "No confounding table was available.",
            ),
        },
        {
            "check_id": "split_perturbation_overlap",
            "question": "Do train and test leak perturbations or target labels across the declared split?",
            "source_section": "split_leakage, combination_leakage",
            "status": max(
                _worst_status(sections.get("split_leakage")),
                _worst_status(sections.get("combination_leakage")),
                key=lambda status: STATUS_RANK.get(status, 0),
            ),
            "message": _message(
                sections.get("split_leakage"),
                "Provide --split to check train/test overlap.",
            ),
        },
        {
            "check_id": "control_balance",
            "question": "Are controls present and balanced across relevant batches, replicates, cell types, or timepoints?",
            "source_section": "control_balance",
            "status": _worst_status(sections.get("control_balance")),
            "message": _message(sections.get("control_balance"), "No control balance table was available."),
        },
        {
            "check_id": "target_effect_detectability",
            "question": "Are claimed target effects detectable in the measured expression matrix?",
            "source_section": "target_effect",
            "status": _worst_status(sections.get("target_effect")),
            "message": _message(sections.get("target_effect"), "No target effect table was available."),
        },
        {
            "check_id": "metadata_only_prediction",
            "question": "Can metadata alone suspiciously predict perturbation identity?",
            "source_section": "metadata_shortcut",
            "status": _worst_status(sections.get("metadata_shortcut")),
            "message": _message(sections.get("metadata_shortcut"), "No metadata-only shortcut table was available."),
        },
        {
            "check_id": "claim_support",
            "question": "Does the supplied split support the claimed unseen perturbation or unseen target generalization?",
            "source_section": "claim_support",
            "status": _worst_status(sections.get("claim_support")),
            "message": _message(
                sections.get("claim_support"),
                "Provide --split and optionally --claim to validate the generalization claim.",
            ),
        },
    ]
    return pd.DataFrame(rows)
