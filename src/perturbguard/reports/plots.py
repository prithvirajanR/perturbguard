from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px


def write_audit_plots(out_dir: str | Path, sections: dict[str, pd.DataFrame]) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    if "cell_count" in sections and {"entity_type", "n_cells"}.issubset(sections["cell_count"].columns):
        table = sections["cell_count"]
        fig = px.histogram(table, x="n_cells", color="entity_type", title="Cell-count support")
        path = out / "cell_count.html"
        fig.write_html(path)
        written["cell_count"] = path
    if "confounding" in sections and {"metadata_variable", "cramers_v"}.issubset(
        sections["confounding"].columns
    ):
        fig = px.bar(
            sections["confounding"],
            x="metadata_variable",
            y="cramers_v",
            color="status",
            title="Perturbation-metadata association",
        )
        path = out / "confounding.html"
        fig.write_html(path)
        written["confounding"] = path
    if "split_balance" in sections:
        fig = px.bar(
            sections["split_balance"],
            x="variable",
            y="max_proportion_difference",
            color="status",
            title="Split balance",
        )
        path = out / "split_balance.html"
        fig.write_html(path)
        written["split_balance"] = path
    return written
