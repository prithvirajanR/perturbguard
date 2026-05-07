from __future__ import annotations

import pandas as pd


STATUS_RANK = {"FAIL": 4, "UNSUPPORTED": 4, "WARNING": 2, "INSUFFICIENT_METADATA": 1, "PASS": 0, "SUPPORTED": 0}


def build_summary(sections: dict[str, pd.DataFrame]) -> dict:
    counts: dict[str, dict[str, int]] = {}
    worst_rank = 0
    worst_status = "PASS"
    for name, table in sections.items():
        if "status" not in table:
            continue
        section_counts = table["status"].value_counts().to_dict()
        counts[name] = {str(k): int(v) for k, v in section_counts.items()}
        for status in section_counts:
            rank = STATUS_RANK.get(str(status), 0)
            if rank > worst_rank:
                worst_rank = rank
                worst_status = str(status)
    return {"overall_status": worst_status, "sections": counts}
