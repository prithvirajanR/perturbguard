from __future__ import annotations

import pandas as pd


DEFAULT_RECOMMENDATIONS = {
    "confounding": "Review batch/donor/replicate randomization before interpreting model performance.",
    "metadata_concentration": "Inspect perturbations concentrated in one batch, donor, or replicate before trusting benchmark results.",
    "control_balance": "Add or restrict to matched controls for affected groups.",
    "cell_count": "Exclude or cautiously interpret low-support perturbations and guides.",
    "metadata_shortcut": "Treat high metadata-only accuracy as evidence of shortcut risk.",
    "target_effect": "Inspect perturbations whose target gene does not move in the expected direction.",
    "guide_consistency": "Review guides with weakly correlated effect vectors.",
    "split_leakage": "Regenerate a split aligned to the intended generalization claim.",
    "claim_support": "Restate the claim or generate a stricter split.",
}


def build_recommendations(sections: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for section, table in sections.items():
        if "status" not in table:
            continue
        flagged = table.loc[table["status"].isin(["FAIL", "WARNING", "UNSUPPORTED"])]
        if flagged.empty:
            continue
        worst = "FAIL" if flagged["status"].isin(["FAIL", "UNSUPPORTED"]).any() else "WARNING"
        rows.append(
            {
                "section": section,
                "status": worst,
                "n_flagged": int(len(flagged)),
                "recommendation": DEFAULT_RECOMMENDATIONS.get(
                    section,
                    "Review flagged rows and decide whether to exclude, repair, or qualify the analysis.",
                ),
            }
        )
    return pd.DataFrame(rows, columns=["section", "status", "n_flagged", "recommendation"])
