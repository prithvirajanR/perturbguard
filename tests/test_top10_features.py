from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from typer.testing import CliRunner

from perturbguard.adversarial.tests import run_adversarial_checks
from perturbguard.benchmark.manifest import check_benchmark_manifest
from perturbguard.benchmark.validation import run_validation_benchmark
from perturbguard.cli.main import app
from perturbguard.design.power import check_power
from perturbguard.evaluate.model import evaluate_predictions
from perturbguard.io.config_wizard import infer_config
from perturbguard.qc.target_mapping import audit_target_mapping
from perturbguard.repair.anndata import repair_anndata
from perturbguard.reports.dataset_card import build_dataset_card
from perturbguard.reports.html import write_html_report
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.streaming.profile import profile_h5ad


runner = CliRunner()


def test_evaluate_predictions_reports_overall_and_group_failures():
    adata = create_synthetic_perturbseq(n_cells=80, random_state=201)
    predictions = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "y_true": adata.obs["perturbation"].astype(str).values,
            "y_pred": adata.obs["perturbation"].astype(str).values,
            "confidence": np.linspace(0.1, 0.9, adata.n_obs),
        }
    )
    predictions.loc[:10, "y_pred"] = "wrong"

    result = evaluate_predictions(adata, predictions, group_columns=["batch"], min_group_size=5)

    assert {"overall", "per_group", "calibration"}.issubset(result)
    assert result["overall"]["n_predictions"].iloc[0] == adata.n_obs
    assert "status" in result["per_group"].columns
    assert result["calibration"]["status"].iloc[0] in {"PASS", "WARNING", "FAIL", "INSUFFICIENT_METADATA"}


def test_evaluate_predictions_reports_perturbation_effect_metrics():
    adata = create_synthetic_perturbseq(n_cells=120, random_state=209)
    controls = adata.obs["is_control"].astype(bool).values
    rows = []
    for perturbation, group in adata.obs.loc[~controls].groupby("perturbation", observed=True):
        idx = adata.obs.index.get_indexer(group.index)
        observed_mean = np.asarray(adata.X[idx]).mean(axis=0)
        predicted_mean = observed_mean.copy()
        row = {
            "perturbation": perturbation,
            "n_predicted_cells": len(group),
        }
        for gene, observed, predicted in zip(adata.var_names[:8], observed_mean[:8], predicted_mean[:8], strict=False):
            row[f"true_{gene}"] = observed
            row[f"pred_{gene}"] = predicted
            row[f"pred_low_{gene}"] = predicted - 0.5
            row[f"pred_high_{gene}"] = predicted + 0.5
        rows.append(row)
    predictions = pd.DataFrame(rows)

    result = evaluate_predictions(
        adata,
        predictions,
        gene_sets={"stress": list(adata.var_names[:4]), "other": list(adata.var_names[4:8])},
        top_k_de=3,
    )

    assert {"effect_recovery", "pathway_recovery", "interval_calibration"}.issubset(result)
    assert result["effect_recovery"]["spearman_effect_corr"].notna().all()
    assert result["effect_recovery"]["top_k_de_recovery"].between(0, 1).all()
    assert result["pathway_recovery"]["pathway_score_correlation"].notna().all()
    assert result["interval_calibration"]["interval_coverage"].between(0, 1).all()


def test_repair_anndata_standardizes_schema_controls_and_duplicate_indices():
    raw = create_synthetic_perturbseq(n_cells=20, random_state=202)
    raw.obs = raw.obs.rename(columns={"perturbation": "condition"})
    raw.obs_names = ["cell"] * raw.n_obs
    raw.var_names = ["gene"] * raw.n_vars

    fixed, actions = repair_anndata(
        raw,
        schema={"perturbation": "condition"},
        control_labels=["control"],
    )

    assert fixed.obs_names.is_unique
    assert fixed.var_names.is_unique
    assert "perturbation" in fixed.obs
    assert "is_control" in fixed.obs
    assert actions["status"].isin(["PASS", "WARNING"]).all()


def test_interactive_report_contains_filter_controls(tmp_path):
    out = tmp_path / "report.html"
    write_html_report(out, {"checks": pd.DataFrame([{"status": "FAIL", "message": "bad"}])})

    html = out.read_text(encoding="utf-8")

    assert "data-status-filter" in html
    assert "Search tables" in html
    assert "summary-card" in html


def test_target_mapping_audits_gene_class_and_configured_drug_targets():
    adata = create_synthetic_perturbseq(n_cells=30, random_state=203)
    adata.obs["target_gene"] = ["TP53", "HDAC inhibitor", "DrugX"] * 10
    mapping = pd.DataFrame(
        [{"perturbation": "DrugX", "target": "MYC", "target_type": "gene", "source": "test"}]
    )

    result = audit_target_mapping(adata, mapping)

    statuses = dict(zip(result["raw_target"], result["target_type"], strict=False))
    assert statuses["TP53"] == "gene"
    assert statuses["HDAC inhibitor"] == "target_class"
    assert statuses["MYC"] == "gene"


def test_power_check_flags_low_cells_and_missing_replicates():
    design = pd.DataFrame(
        {
            "batch": ["b1", "b1"],
            "perturbation": ["control", "TP53"],
            "replicate": ["r1", "r1"],
            "n_cells_planned": [100, 12],
        }
    )

    result = check_power(design, min_cells_per_perturbation=50, min_replicates=2)

    assert result["status"].isin(["FAIL", "WARNING"]).any()
    assert {"planned_cells_per_perturbation", "replicate_support"}.issubset(set(result["check_id"]))


def test_benchmark_manifest_validates_paths_and_claim(tmp_path):
    adata = create_synthetic_perturbseq(n_cells=60, random_state=204)
    data = tmp_path / "data.h5ad"
    split = tmp_path / "split.csv"
    manifest = tmp_path / "benchmark.yaml"
    adata.write_h5ad(data)
    pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": "train"}).to_csv(split, index=False)
    manifest.write_text(
        yaml.safe_dump(
            {
                "dataset": str(data),
                "split": str(split),
                "claim": "unseen_perturbation",
                "model": {"name": "demo"},
                "metrics": ["accuracy"],
            }
        ),
        encoding="utf-8",
    )

    result = check_benchmark_manifest(manifest)

    assert result.loc[result["check_id"].eq("manifest_file"), "status"].iloc[0] == "PASS"
    assert result.loc[result["check_id"].eq("claim_support"), "status"].iloc[0] in {
        "UNSUPPORTED",
        "INSUFFICIENT_METADATA",
    }


def test_validation_benchmark_scores_expected_audit_findings(tmp_path):
    adata = create_synthetic_perturbseq(scenario="batch_confounded", n_cells=180, random_state=211)
    data = tmp_path / "batch_confounded.h5ad"
    expected = tmp_path / "expected_findings.csv"
    manifest = tmp_path / "validation.yaml"
    adata.write_h5ad(data)
    pd.DataFrame(
        [
            {
                "case": "batch_confounded",
                "section": "metadata_concentration",
                "match_column": "metadata_variable",
                "match_value": "batch",
                "expected_status": "FAIL",
            },
            {
                "case": "batch_confounded",
                "section": "metadata_shortcut",
                "expected_status": "FAIL",
            },
        ]
    ).to_csv(expected, index=False)
    manifest.write_text(
        yaml.safe_dump(
            {
                "cases": [{"name": "batch_confounded", "dataset": str(data)}],
                "expected_findings": str(expected),
            }
        ),
        encoding="utf-8",
    )

    result = run_validation_benchmark(manifest)

    assert {"case", "section", "expected_status", "observed_status", "status"}.issubset(result.columns)
    assert result["status"].eq("PASS").all()


def test_cli_validation_benchmark_writes_expected_finding_results(tmp_path):
    adata = create_synthetic_perturbseq(scenario="batch_confounded", n_cells=180, random_state=212)
    data = tmp_path / "data.h5ad"
    expected = tmp_path / "expected.csv"
    manifest = tmp_path / "validation.yaml"
    out = tmp_path / "validation"
    adata.write_h5ad(data)
    pd.DataFrame(
        [
            {
                "case": "case1",
                "section": "metadata_concentration",
                "match_column": "metadata_variable",
                "match_value": "batch",
                "expected_status": "FAIL",
            }
        ]
    ).to_csv(expected, index=False)
    manifest.write_text(
        yaml.safe_dump(
            {
                "cases": [{"name": "case1", "dataset": str(data)}],
                "expected_findings": str(expected),
            }
        ),
        encoding="utf-8",
    )

    result = runner.invoke(app, ["expected-findings", "--manifest", str(manifest), "--out", str(out)])

    assert result.exit_code == 0, result.stdout
    assert (out / "validation_results.csv").exists()
    assert pd.read_csv(out / "validation_results.csv")["status"].eq("PASS").all()


def test_profile_h5ad_reports_backed_shape_without_crashing(tmp_path):
    adata = create_synthetic_perturbseq(n_cells=25, random_state=205)
    data = tmp_path / "profile.h5ad"
    adata.write_h5ad(data)

    result = profile_h5ad(data)

    assert result.loc[result["check_id"].eq("backed_load"), "status"].iloc[0] == "PASS"
    assert result.loc[result["check_id"].eq("n_cells"), "value"].iloc[0] == 25


def test_adversarial_checks_detect_metadata_shortcuts():
    adata = create_synthetic_perturbseq(scenario="batch_confounded", n_cells=100, random_state=206)

    result = run_adversarial_checks(adata, min_cells_per_class=3, n_splits=3)

    assert "metadata_predicts_perturbation" in set(result["check_id"])
    assert result["status"].isin(["WARNING", "FAIL"]).any()


def test_infer_config_finds_common_column_aliases():
    obs = pd.DataFrame(
        {
            "condition": ["control", "TP53"],
            "batch_id": ["b1", "b2"],
            "guide": ["ntc", "g1"],
        }
    )
    adata = ad.AnnData(np.ones((2, 2)), obs=obs, var=pd.DataFrame(index=["TP53", "MYC"]))

    config = infer_config(adata)

    assert config["schema"]["perturbation"] == "condition"
    assert config["schema"]["batch"] == "batch_id"
    assert config["schema"]["guide_id"] == "guide"


def test_dataset_card_contains_claims_and_limitations():
    adata = create_synthetic_perturbseq(n_cells=30, random_state=207)
    sections = {"dataset_validation": pd.DataFrame([{"status": "PASS", "message": "ok"}])}

    card = build_dataset_card(adata, sections, dataset_name="demo")

    assert "# Dataset Card: demo" in card
    assert "Known Limitations" in card
    assert "n_cells" in card


def test_cli_exposes_new_commands(tmp_path):
    adata = create_synthetic_perturbseq(n_cells=40, random_state=208)
    data = tmp_path / "data.h5ad"
    pred = tmp_path / "pred.csv"
    adata.write_h5ad(data)
    pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "y_true": adata.obs["perturbation"].astype(str).values,
            "y_pred": adata.obs["perturbation"].astype(str).values,
        }
    ).to_csv(pred, index=False)

    commands = runner.invoke(app, ["--help"])
    assert commands.exit_code == 0
    for command in [
        "evaluate",
        "repair",
        "infer-config",
        "target-map",
        "power-check",
        "benchmark-check",
        "expected-findings",
        "validation-benchmark",
        "profile-large",
        "adversarial-check",
        "dataset-card",
    ]:
        assert command in commands.stdout

    result = runner.invoke(app, ["evaluate", "--data", str(data), "--predictions", str(pred), "--out", str(tmp_path / "eval")])
    assert result.exit_code == 0, result.stdout
    assert json.loads((tmp_path / "eval" / "summary.json").read_text(encoding="utf-8"))["overall_status"] in {
        "PASS",
        "WARNING",
        "FAIL",
    }


def test_cli_evaluate_writes_effect_recovery_tables_for_expression_predictions(tmp_path):
    adata = create_synthetic_perturbseq(n_cells=80, random_state=210)
    data = tmp_path / "data.h5ad"
    pred = tmp_path / "effect_predictions.csv"
    gene_sets = tmp_path / "gene_sets.yaml"
    adata.write_h5ad(data)
    rows = []
    for perturbation, group in adata.obs.loc[~adata.obs["is_control"].astype(bool)].groupby("perturbation", observed=True):
        idx = adata.obs.index.get_indexer(group.index)
        observed_mean = np.asarray(adata.X[idx]).mean(axis=0)
        row = {"perturbation": perturbation}
        for gene, observed in zip(adata.var_names[:6], observed_mean[:6], strict=False):
            row[f"true_{gene}"] = observed
            row[f"pred_{gene}"] = observed
        rows.append(row)
    pd.DataFrame(rows).to_csv(pred, index=False)
    gene_sets.write_text(yaml.safe_dump({"set_a": list(adata.var_names[:3]), "set_b": list(adata.var_names[3:6])}), encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "evaluate",
            "--data",
            str(data),
            "--predictions",
            str(pred),
            "--gene-sets",
            str(gene_sets),
            "--out",
            str(tmp_path / "eval"),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert (tmp_path / "eval" / "tables" / "effect_recovery.csv").exists()
    assert (tmp_path / "eval" / "tables" / "pathway_recovery.csv").exists()
