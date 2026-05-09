from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


STATUS_RANK = {
    "PASS": 0,
    "SUPPORTED": 0,
    "INSUFFICIENT_METADATA": 1,
    "WARNING": 2,
    "FAIL": 3,
    "UNSUPPORTED": 3,
}


def _resolve(base: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else base / path


def _worst_status(table: pd.DataFrame) -> str:
    if "status" not in table or table.empty:
        return "INSUFFICIENT_METADATA"
    statuses = table["status"].dropna().astype(str)
    if statuses.empty:
        return "INSUFFICIENT_METADATA"
    return max(statuses, key=lambda status: STATUS_RANK.get(status, -1))


def _observed_status(sections: dict[str, pd.DataFrame], expected: pd.Series) -> str:
    section = str(expected["section"])
    if section not in sections:
        return "MISSING_SECTION"
    table = sections[section]
    match_column = expected.get("match_column")
    match_value = expected.get("match_value")
    if pd.notna(match_column) and pd.notna(match_value) and str(match_column) in table:
        table = table.loc[table[str(match_column)].astype(str).eq(str(match_value))]
    return _worst_status(table)


def run_validation_benchmark(manifest_path: str | Path) -> pd.DataFrame:
    from perturbguard.cli.main import _audit_sections

    manifest_path = Path(manifest_path)
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    base = manifest_path.parent
    expected_path = _resolve(base, manifest["expected_findings"])
    expected = pd.read_csv(expected_path)
    sections_by_case = {}
    for case in manifest.get("cases", []):
        name = str(case["name"])
        data_path = _resolve(base, case["dataset"])
        config_path = _resolve(base, case["config"]) if case.get("config") else None
        split_path = _resolve(base, case["split"]) if case.get("split") else None
        sections_by_case[name] = _audit_sections(
            data_path,
            config=config_path,
            split=split_path,
            claim_name=case.get("claim"),
        )
    rows = []
    for _, row in expected.iterrows():
        case = str(row["case"])
        observed = (
            _observed_status(sections_by_case[case], row)
            if case in sections_by_case
            else "MISSING_CASE"
        )
        expected_status = str(row["expected_status"])
        rows.append(
            {
                "case": case,
                "section": str(row["section"]),
                "match_column": "" if pd.isna(row.get("match_column")) else str(row.get("match_column")),
                "match_value": "" if pd.isna(row.get("match_value")) else str(row.get("match_value")),
                "expected_status": expected_status,
                "observed_status": observed,
                "status": "PASS" if observed == expected_status else "FAIL",
                "message": f"Expected {expected_status}, observed {observed}.",
            }
        )
    return pd.DataFrame(rows)
