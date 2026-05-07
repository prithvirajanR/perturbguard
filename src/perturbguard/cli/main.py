from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd
import typer
import yaml

from perturbguard.adversarial.tests import run_adversarial_checks
from perturbguard.benchmark.manifest import check_benchmark_manifest
from perturbguard.claims.claim_checker import check_claim
from perturbguard.compare.datasets import compare_datasets
from perturbguard.design.checker import check_design
from perturbguard.design.power import check_power
from perturbguard.evaluate.model import evaluate_predictions
from perturbguard.io.config import load_config
from perturbguard.io.config_wizard import infer_config
from perturbguard.io.schema import apply_schema
from perturbguard.io.split_loader import load_split
from perturbguard.leakage.combination_leakage import check_combination_leakage
from perturbguard.leakage.split_leakage import check_split_leakage
from perturbguard.leakage.split_balance import evaluate_split_balance
from perturbguard.qc.confounding import detect_confounding
from perturbguard.qc.cell_count import check_cell_counts
from perturbguard.qc.control_balance import check_control_balance
from perturbguard.qc.dataset_validator import validate_anndata, validate_h5ad_file
from perturbguard.qc.guide_consistency import check_guide_consistency
from perturbguard.qc.metadata_shortcut import metadata_shortcut_baseline
from perturbguard.qc.perturbation_summary import summarize_perturbations
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.qc.target_mapping import audit_target_mapping
from perturbguard.repair.anndata import repair_anndata
from perturbguard.reports.dataset_card import build_dataset_card
from perturbguard.reports.html import write_html_report
from perturbguard.reports.plots import write_audit_plots
from perturbguard.reports.recommendations import build_recommendations
from perturbguard.reports.summary import build_summary
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split
from perturbguard.streaming.profile import profile_h5ad

app = typer.Typer(help="Audit perturbation datasets, splits, and claims.")


@app.command()
def simulate(scenario: str = "clean", out: Path = typer.Option(...), n_cells: int = 360, n_genes: int = 60):
    out.parent.mkdir(parents=True, exist_ok=True)
    create_synthetic_perturbseq(scenario=scenario, n_cells=n_cells, n_genes=n_genes).write_h5ad(out)
    typer.echo(f"Wrote {out}")


@app.command()
def validate(data: Path = typer.Option(...), out: Path | None = None, config: Path | None = None):
    cfg = load_config(config)
    if cfg.schema:
        try:
            df = validate_anndata(apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels")))
        except Exception:
            df = validate_h5ad_file(data)
    else:
        df = validate_h5ad_file(data)
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
    typer.echo(df.to_string(index=False))


def _audit_sections(
    data: Path,
    config: Path | None = None,
    split: Path | None = None,
    claim_name: str | None = None,
):
    cfg = load_config(config)
    adata = apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels"))
    dataset_validation = validate_anndata(adata)
    if dataset_validation.loc[
        dataset_validation["check_id"].eq("obs_perturbation"), "status"
    ].eq("FAIL").any():
        return {
            "dataset_validation": dataset_validation,
            "audit_status": pd.DataFrame(
                [
                    {
                        "status": "FAIL",
                        "message": "Audit stopped because required perturbation metadata is missing.",
                    }
                ]
            ),
        }
    sections = {
        "dataset_validation": dataset_validation,
        "perturbation_summary": summarize_perturbations(adata),
        "control_balance": check_control_balance(adata),
        "cell_count": check_cell_counts(adata),
        "confounding": detect_confounding(adata),
        "metadata_shortcut": metadata_shortcut_baseline(adata),
        "target_effect": check_target_effect(adata),
        "target_mapping": audit_target_mapping(adata),
        "guide_consistency": check_guide_consistency(adata),
    }
    if split is not None:
        split_df = load_split(split)
        sections["split_leakage"] = check_split_leakage(
            adata,
            split_df,
            claim=claim_name or "unseen_perturbation",
        )
        sections["combination_leakage"] = check_combination_leakage(adata, split_df)
        sections["split_balance"] = evaluate_split_balance(adata, split_df)
        if claim_name is not None:
            sections["claim_support"] = check_claim(adata, split_df, claim_name)
    return sections


@app.command()
def audit(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
    split: Path | None = None,
    claim_name: str | None = typer.Option(None, "--claim"),
    fail_on: str = typer.Option("never", help="Exit nonzero on audit status: never, fail, warning."),
):
    out.mkdir(parents=True, exist_ok=True)
    tables = out / "tables"
    tables.mkdir(exist_ok=True)
    sections = _audit_sections(data, config, split, claim_name)
    recommendations = build_recommendations(sections)
    sections["recommendations"] = recommendations
    for name, df in sections.items():
        df.to_csv(tables / f"{name}.csv", index=False)
    plot_paths = write_audit_plots(out / "plots", sections)
    write_html_report(out / "report.html", sections, plot_paths)
    summary = build_summary(sections)
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(f"Wrote audit report to {out}")
    normalized_fail_on = fail_on.lower()
    if normalized_fail_on not in {"never", "fail", "warning"}:
        raise typer.BadParameter("--fail-on must be one of: never, fail, warning")
    if normalized_fail_on == "fail" and summary["overall_status"] in {"FAIL", "UNSUPPORTED"}:
        raise typer.Exit(1)
    if normalized_fail_on == "warning" and summary["overall_status"] in {"FAIL", "UNSUPPORTED", "WARNING"}:
        raise typer.Exit(1)


@app.command(name="split")
def split_cmd(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    strategy: str = "random",
    test_size: float = 0.2,
    val_size: float = 0.0,
    random_state: int = 0,
    metadata_column: str | None = None,
    config: Path | None = None,
):
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config)
    split_df = generate_split(
        apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels")),
        strategy=strategy,
        test_size=test_size,
        val_size=val_size,
        random_state=random_state,
        metadata_column=metadata_column,
    )
    split_df.to_csv(out / "split.csv", index=False)
    typer.echo(f"Wrote {out / 'split.csv'}")


@app.command()
def claim(
    data: Path = typer.Option(...),
    split: Path = typer.Option(...),
    claim_name: str = typer.Option("unseen_perturbation", "--claim"),
    out: Path = typer.Option(...),
    config: Path | None = None,
):
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config)
    result = check_claim(
        apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels")),
        load_split(split),
        claim_name,
    )
    result.to_csv(out / "claim_support.csv", index=False)
    typer.echo(result.to_string(index=False))


@app.command()
def evaluate(
    data: Path = typer.Option(...),
    predictions: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
    group_column: list[str] = typer.Option(None, "--group-column"),
):
    out.mkdir(parents=True, exist_ok=True)
    tables = out / "tables"
    tables.mkdir(exist_ok=True)
    cfg = load_config(config)
    adata = apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels"))
    sections = evaluate_predictions(adata, pd.read_csv(predictions), group_columns=group_column or None)
    for name, df in sections.items():
        df.to_csv(tables / f"{name}.csv", index=False)
    write_html_report(out / "report.html", sections)
    summary = build_summary(sections)
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(f"Wrote model evaluation to {out}")


@app.command()
def repair(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
):
    out.parent.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config)
    fixed, actions = repair_anndata(
        ad.read_h5ad(data),
        schema=cfg.schema,
        control_labels=cfg.controls.get("labels"),
    )
    fixed.write_h5ad(out)
    actions_path = out.with_suffix(".repair_actions.csv")
    actions.to_csv(actions_path, index=False)
    typer.echo(f"Wrote repaired AnnData to {out}")
    typer.echo(f"Wrote repair actions to {actions_path}")


@app.command(name="infer-config")
def infer_config_cmd(data: Path = typer.Option(...), out: Path = typer.Option(...)):
    out.parent.mkdir(parents=True, exist_ok=True)
    config = infer_config(ad.read_h5ad(data, backed="r"))
    out.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    typer.echo(f"Wrote inferred config to {out}")


@app.command(name="target-map")
def target_map_cmd(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
    mapping: Path | None = None,
):
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config)
    target_map = pd.read_csv(mapping) if mapping else None
    result = audit_target_mapping(
        apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels")),
        target_map,
    )
    result.to_csv(out / "target_mapping.csv", index=False)
    write_html_report(out / "report.html", {"target_mapping": result})
    typer.echo(f"Wrote target mapping audit to {out}")


@app.command(name="compare-datasets")
def compare_datasets_cmd(data: list[Path] = typer.Option(...), out: Path = typer.Option(...)):
    out.mkdir(parents=True, exist_ok=True)
    datasets = {path.stem: ad.read_h5ad(path) for path in data}
    result = compare_datasets(datasets)
    result.to_csv(out / "dataset_comparison.csv", index=False)
    write_html_report(out / "report.html", {"dataset_comparison": result})
    typer.echo(f"Wrote dataset comparison to {out}")


@app.command(name="design-check")
def design_check_cmd(design: Path = typer.Option(...), out: Path = typer.Option(...)):
    out.mkdir(parents=True, exist_ok=True)
    if design.suffix.lower() in {".yaml", ".yml"}:
        raw = load_config(design).raw
        table = pd.DataFrame(raw.get("design", raw if isinstance(raw, list) else []))
    else:
        table = pd.read_csv(design, sep="\t" if design.suffix.lower() == ".tsv" else ",")
    result = check_design(table)
    result.to_csv(out / "design_validation.csv", index=False)
    write_html_report(out / "report.html", {"design_validation": result})
    typer.echo(f"Wrote design report to {out}")


@app.command(name="power-check")
def power_check_cmd(design: Path = typer.Option(...), out: Path = typer.Option(...)):
    out.mkdir(parents=True, exist_ok=True)
    if design.suffix.lower() in {".yaml", ".yml"}:
        raw = load_config(design).raw
        table = pd.DataFrame(raw.get("design", raw if isinstance(raw, list) else []))
    else:
        table = pd.read_csv(design, sep="\t" if design.suffix.lower() == ".tsv" else ",")
    result = check_power(table)
    result.to_csv(out / "power_validation.csv", index=False)
    write_html_report(out / "report.html", {"power_validation": result})
    typer.echo(f"Wrote power report to {out}")


@app.command(name="benchmark-check")
def benchmark_check_cmd(manifest: Path = typer.Option(...), out: Path = typer.Option(...)):
    out.mkdir(parents=True, exist_ok=True)
    result = check_benchmark_manifest(manifest)
    result.to_csv(out / "benchmark_validation.csv", index=False)
    write_html_report(out / "report.html", {"benchmark_validation": result})
    typer.echo(f"Wrote benchmark validation to {out}")


@app.command(name="profile-large")
def profile_large_cmd(data: Path = typer.Option(...), out: Path = typer.Option(...)):
    out.mkdir(parents=True, exist_ok=True)
    result = profile_h5ad(data)
    result.to_csv(out / "large_profile.csv", index=False)
    write_html_report(out / "report.html", {"large_profile": result})
    typer.echo(f"Wrote large-file profile to {out}")


@app.command(name="adversarial-check")
def adversarial_check_cmd(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
    split: Path | None = None,
):
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config)
    split_df = load_split(split) if split else None
    result = run_adversarial_checks(
        apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels")),
        split_df,
    )
    result.to_csv(out / "adversarial_checks.csv", index=False)
    write_html_report(out / "report.html", {"adversarial_checks": result})
    typer.echo(f"Wrote adversarial checks to {out}")


@app.command(name="dataset-card")
def dataset_card_cmd(
    data: Path = typer.Option(...),
    out: Path = typer.Option(...),
    config: Path | None = None,
    name: str | None = None,
):
    out.parent.mkdir(parents=True, exist_ok=True)
    sections = _audit_sections(data, config)
    cfg = load_config(config)
    adata = apply_schema(ad.read_h5ad(data), cfg.schema, cfg.controls.get("labels"))
    card = build_dataset_card(adata, sections, dataset_name=name or data.stem)
    out.write_text(card, encoding="utf-8")
    typer.echo(f"Wrote dataset card to {out}")
