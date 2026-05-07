from __future__ import annotations

import pandas as pd


def result(check_id: str, status: str, message: str, affected_column: str = "", suggested_fix: str = "") -> dict[str, str]:
    return {
        "check_id": check_id,
        "status": status,
        "message": message,
        "affected_column": affected_column,
        "suggested_fix": suggested_fix,
    }


def frame(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)
