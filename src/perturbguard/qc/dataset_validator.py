from __future__ import annotations

import pandas as pd
from anndata import AnnData
import anndata as ad
from pathlib import Path

from perturbguard.io.schema import Schema
from perturbguard.utils.status import frame, result


def validate_anndata(adata: AnnData, schema: Schema | None = None) -> pd.DataFrame:
    schema = schema or Schema()
    rows = []
    rows.append(result("matrix_exists", "PASS" if adata.X is not None else "FAIL", "Expression matrix is present." if adata.X is not None else "Expression matrix is missing.", "X", "Provide an AnnData object with X."))
    rows.append(result("n_cells", "PASS" if adata.n_obs > 0 else "FAIL", f"Dataset has {adata.n_obs} cells.", "obs", "Provide at least one cell."))
    rows.append(result("n_genes", "PASS" if adata.n_vars > 0 else "FAIL", f"Dataset has {adata.n_vars} genes.", "var", "Provide at least one gene."))
    rows.append(result("obs_columns_unique", "PASS" if adata.obs.columns.is_unique else "FAIL", "Observation metadata columns are unique." if adata.obs.columns.is_unique else "Observation metadata column names are duplicated.", "obs", "Make obs column names unique."))
    if not adata.obs.columns.is_unique:
        return frame(rows)
    pert_col = schema.resolve_column(adata.obs, "perturbation")
    rows.append(result("obs_perturbation", "PASS" if pert_col else "FAIL", "Perturbation column is present." if pert_col else "Perturbation column is missing.", schema.mapping.get("perturbation", "perturbation"), "Map or add a perturbation column."))
    if pert_col:
        try:
            has_controls = bool(schema.get_control_mask(adata.obs).any())
        except KeyError:
            has_controls = False
        rows.append(result("controls_exist", "PASS" if has_controls else "FAIL", "Control cells are present." if has_controls else "No controls found.", "is_control", "Add is_control or configure control labels."))
    rows.append(result("var_names_unique", "PASS" if adata.var_names.is_unique else "WARNING", "Gene names are unique." if adata.var_names.is_unique else "Gene names are duplicated.", "var_names", "Make gene names unique."))
    rows.append(result("obs_names_unique", "PASS" if adata.obs_names.is_unique else "WARNING", "Cell names are unique." if adata.obs_names.is_unique else "Cell names are duplicated.", "obs_names", "Make cell names unique before split alignment."))
    for col in ["batch", "replicate"]:
        actual = schema.resolve_column(adata.obs, col)
        rows.append(result(f"obs_{col}", "PASS" if actual else "INSUFFICIENT_METADATA", f"{col} column is present." if actual else f"{col} column is unavailable.", col, f"Add or map {col} metadata."))
    for col in ["dose", "timepoint"]:
        actual = schema.resolve_column(adata.obs, col)
        if actual:
            mixed = adata.obs[actual].map(lambda x: type(x).__name__).nunique(dropna=True) > 1
            rows.append(result(f"metadata_type_{col}", "WARNING" if mixed else "PASS", f"{col} metadata has mixed Python value types." if mixed else f"{col} metadata type is consistent.", actual, "Normalize metadata values to a single dtype."))
    return frame(rows)


def validate_h5ad_file(data: str | Path | AnnData, schema: Schema | None = None) -> pd.DataFrame:
    if isinstance(data, AnnData):
        rows = [result("file_load", "PASS", "AnnData object is already loaded.", "", "")]
        return pd.concat([pd.DataFrame(rows), validate_anndata(data, schema)], ignore_index=True)
    try:
        adata = ad.read_h5ad(data)
    except Exception as exc:
        return pd.DataFrame(
            [
                result(
                    "file_load",
                    "FAIL",
                    f"Failed to load AnnData file: {exc}",
                    str(data),
                    "Provide a valid .h5ad file.",
                )
            ]
        )
    rows = [result("file_load", "PASS", "AnnData file loads successfully.", str(data), "")]
    return pd.concat([pd.DataFrame(rows), validate_anndata(adata, schema)], ignore_index=True)
