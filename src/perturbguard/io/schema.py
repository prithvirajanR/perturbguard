from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd
from anndata import AnnData
from pandas.api.types import is_bool_dtype
from perturbguard.utils.controls import control_mask


DEFAULT_CONTROL_LABELS = {"control", "ctrl", "ntc", "non-targeting", "dmso", "vehicle"}


@dataclass
class Schema:
    mapping: dict[str, str] = field(default_factory=dict)

    def resolve_column(self, df: pd.DataFrame, canonical: str) -> str | None:
        candidate = self.mapping.get(canonical, canonical)
        return candidate if candidate in df.columns else None

    def require(self, df: pd.DataFrame, canonical: str) -> str:
        column = self.resolve_column(df, canonical)
        if column is None:
            raise KeyError(f"Missing required column for {canonical!r}")
        return column

    def optional(self, df: pd.DataFrame, canonical: str) -> str | None:
        return self.resolve_column(df, canonical)

    def get_control_mask(self, df: pd.DataFrame, control_labels: list[str] | None = None) -> pd.Series:
        control_col = self.resolve_column(df, "is_control")
        if control_col:
            return control_mask(df.rename(columns={control_col: "is_control"}))
        perturbation_col = self.require(df, "perturbation")
        labels = {x.lower() for x in (control_labels or sorted(DEFAULT_CONTROL_LABELS))}
        return df[perturbation_col].astype(str).str.lower().isin(labels)


def apply_schema(
    adata: AnnData,
    mapping: dict[str, str] | None = None,
    control_labels: list[str] | None = None,
) -> AnnData:
    mapping = mapping or {}
    if not mapping:
        return adata.to_memory().copy() if getattr(adata, "isbacked", False) else adata
    mapped = adata.to_memory().copy() if getattr(adata, "isbacked", False) else adata.copy()
    for canonical, source in mapping.items():
        if source in mapped.obs and canonical not in mapped.obs:
            mapped.obs[canonical] = mapped.obs[source]
        elif source in mapped.obs and canonical in mapped.obs and canonical != source:
            mapped.obs[canonical] = mapped.obs[source]
    if "is_control" in mapping and "is_control" in mapped.obs and control_labels:
        if not is_bool_dtype(mapped.obs["is_control"]):
            labels = {str(label).strip().lower() for label in control_labels}
            mapped.obs["is_control"] = mapped.obs["is_control"].astype(str).str.strip().str.lower().isin(labels)
    return mapped
