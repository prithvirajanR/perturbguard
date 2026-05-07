from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from perturbguard.claims.claim_checker import check_claim
from perturbguard.cli.main import app
from perturbguard.io.schema import Schema
from perturbguard.leakage.combination_leakage import categorize_combination_leakage
from perturbguard.leakage.split_leakage import check_split_leakage
from perturbguard.qc.confounding import detect_confounding
from perturbguard.qc.control_balance import check_control_balance
from perturbguard.qc.dataset_validator import validate_anndata
from perturbguard.qc.guide_consistency import check_guide_consistency
from perturbguard.qc.metadata_shortcut import metadata_shortcut_baseline
from perturbguard.qc.perturbation_summary import summarize_perturbations
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.reports.html import write_html_report
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split


runner = CliRunner()


def test_cli_help_runs():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "audit" in result.stdout
    assert "simulate" in result.stdout


def test_synthetic_scenarios_have_expected_metadata():
    adata = create_synthetic_perturbseq(scenario="clean", n_cells=240, n_genes=40, random_state=7)
    assert adata.X.shape == (240, 40)
    for column in [
        "perturbation",
        "is_control",
        "target_gene",
        "guide_id",
        "perturbation_type",
        "batch",
        "replicate",
        "cell_type",
    ]:
        assert column in adata.obs
    assert adata.obs["is_control"].any()


def test_schema_resolves_mapped_columns_and_infers_controls():
    obs = pd.DataFrame({"condition": ["ctrl", "TP53"], "batch_id": ["b1", "b1"]})
    schema = Schema({"perturbation": "condition", "batch": "batch_id"})
    assert schema.require(obs, "perturbation") == "condition"
    controls = schema.get_control_mask(obs, control_labels=["ctrl"])
    assert controls.tolist() == [True, False]


def test_validator_returns_structured_fail_for_missing_perturbation():
    adata = create_synthetic_perturbseq()
    adata.obs = adata.obs.drop(columns=["perturbation"])
    results = validate_anndata(adata)
    row = results.loc[results["check_id"] == "obs_perturbation"].iloc[0]
    assert row["status"] == "FAIL"
    assert set(["check_id", "status", "message", "affected_column", "suggested_fix"]).issubset(results.columns)


def test_qc_modules_flag_synthetic_problems():
    confounded = create_synthetic_perturbseq(scenario="batch_confounded", random_state=1)
    confounding = detect_confounding(confounded)
    assert confounding.loc[confounding["metadata_variable"] == "batch", "status"].iloc[0] == "FAIL"

    control_balance = check_control_balance(confounded)
    assert "status" in control_balance.columns

    shortcut = metadata_shortcut_baseline(confounded, n_splits=3, random_state=1)
    assert shortcut["status"].iloc[0] in {"WARNING", "FAIL"}

    low_count = create_synthetic_perturbseq(scenario="low_cell_count", random_state=2)
    summary = summarize_perturbations(low_count, min_cells_warning=50, min_cells_fail=20)
    assert (summary["status"] != "PASS").any()


def test_target_effect_and_guide_consistency_detect_problem_scenarios():
    failed = create_synthetic_perturbseq(scenario="failed_knockdown", random_state=3)
    target = check_target_effect(failed)
    assert (target["status"] == "FAIL").any()

    inconsistent = create_synthetic_perturbseq(scenario="guide_inconsistent", random_state=4)
    guide = check_guide_consistency(inconsistent)
    assert (guide["status"] != "PASS").any()


def test_split_leakage_combination_claim_and_generator():
    adata = create_synthetic_perturbseq(scenario="leaky_combo", random_state=5)
    split = generate_split(adata, strategy="leave_target_gene_out", test_size=0.25, random_state=5)
    leakage = check_split_leakage(adata, split, claim="unseen_target_gene")
    assert leakage.loc[leakage["leakage_type"] == "target_gene_overlap", "status"].iloc[0] == "PASS"

    combo = categorize_combination_leakage(
        train_perturbations=["A", "B", "A+B"],
        test_perturbations=["A+B", "A+C", "D+E"],
    )
    assert set(combo["category"]) >= {
        "EXACT_LEAKAGE",
        "PARTIAL_COMPONENT_HOLDOUT",
        "STRICT_UNSEEN_COMBINATION",
    }

    claim = check_claim(adata, split, claim="unseen_target_gene")
    assert claim["status"].iloc[0] in {"SUPPORTED", "INSUFFICIENT_METADATA"}


def test_cli_simulate_audit_claim_and_report(tmp_path: Path):
    data_path = tmp_path / "demo.h5ad"
    out_dir = tmp_path / "audit"

    result = runner.invoke(app, ["simulate", "--scenario", "clean", "--out", str(data_path)])
    assert result.exit_code == 0, result.stdout
    assert data_path.exists()

    result = runner.invoke(app, ["audit", "--data", str(data_path), "--out", str(out_dir)])
    assert result.exit_code == 0, result.stdout
    assert (out_dir / "report.html").exists()
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "tables" / "dataset_validation.csv").exists()

    split_dir = tmp_path / "split"
    result = runner.invoke(
        app,
        [
            "split",
            "--data",
            str(data_path),
            "--strategy",
            "leave-target-gene-out",
            "--out",
            str(split_dir),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert (split_dir / "split.csv").exists()

    adata = create_synthetic_perturbseq()
    report = write_html_report(out_dir / "manual_report.html", {"dataset_validation": validate_anndata(adata)})
    assert report.exists()
    assert "Dataset Validation" in report.read_text(encoding="utf-8")
