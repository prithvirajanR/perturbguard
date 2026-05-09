from pathlib import Path
import tomllib
from typer.testing import CliRunner

import perturbguard
from perturbguard.cli.main import app
from perturbguard.leakage.split_balance import evaluate_split_balance
from perturbguard.qc.dataset_validator import validate_h5ad_file
from perturbguard.qc.target_mapping import classify_target
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.reports.plots import write_audit_plots
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split


runner = CliRunner()


def test_validate_h5ad_file_returns_fail_for_corrupted_file(tmp_path: Path):
    bad = tmp_path / "bad.h5ad"
    bad.write_bytes(b"not an h5ad")

    result = validate_h5ad_file(bad)

    assert result.loc[result["check_id"].eq("file_load"), "status"].iloc[0] == "FAIL"


def test_validator_reports_duplicate_obs_names_and_mixed_metadata():
    adata = create_synthetic_perturbseq(random_state=81)
    adata.obs_names = ["dup"] * adata.n_obs
    adata.obs["dose"] = ["1", 2, None] * (adata.n_obs // 3)

    result = validate_h5ad_file(adata)

    assert result.loc[result["check_id"].eq("obs_names_unique"), "status"].iloc[0] == "WARNING"
    assert "metadata_type_dose" in set(result["check_id"])


def test_target_mapping_distinguishes_gene_pathway_drug_and_missing():
    assert classify_target("TP53", {"TP53", "MYC"}).target_type == "gene"
    assert classify_target("JAK/STAT", {"TP53", "MYC"}).target_type == "pathway_or_class"
    assert classify_target("HDAC", {"TP53", "MYC"}).target_type == "target_class"
    assert classify_target("", {"TP53", "MYC"}).target_type == "missing"


def test_target_effect_reports_target_type_for_unmapped_drug_or_variant_targets():
    adata = create_synthetic_perturbseq(random_state=811)
    adata.obs.loc[~adata.obs["is_control"].astype(bool), "target_gene"] = "HDAC"

    results = check_target_effect(adata)

    assert "target_type" in results.columns
    assert results["target_type"].eq("target_class").any()


def test_split_balance_report_and_plots_are_written(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=82)
    split = generate_split(adata, strategy="random", random_state=82)

    balance = evaluate_split_balance(adata, split, variables=["batch", "cell_type"])
    plots = write_audit_plots(tmp_path, {"split_balance": balance})

    assert {"variable", "max_proportion_difference", "status"}.issubset(balance.columns)
    assert (tmp_path / "split_balance.html").exists()
    assert plots["split_balance"].exists()


def test_audit_writes_plots_directory(tmp_path: Path):
    data = tmp_path / "data.h5ad"
    out = tmp_path / "audit"
    create_synthetic_perturbseq(random_state=83).write_h5ad(data)

    result = runner.invoke(app, ["audit", "--data", str(data), "--out", str(out)])

    assert result.exit_code == 0, result.stdout
    assert (out / "plots").exists()
    assert any((out / "plots").glob("*.html"))


def test_production_scaffolding_files_exist_and_metadata_is_versioned():
    assert Path(".github/workflows/ci.yml").exists()
    assert Path("CHANGELOG.md").exists()
    assert Path("LICENSE").exists()
    assert Path("docs/assumptions.md").exists()
    assert Path("docs/input_contract.md").exists()


def test_package_version_and_maturity_metadata_are_consistent():
    pyproject = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]

    assert perturbguard.__version__ == project["version"]
    assert "Development Status :: 4 - Beta" in project["classifiers"]
    assert "Development Status :: 5 - Production/Stable" not in project["classifiers"]

    readme = Path("README.md").read_text(encoding="utf-8")
    assert f"Current source version: `{project['version']}`" in readme
    assert "Current stable release" not in readme


def test_docs_do_not_present_smoke_tests_as_scientific_validation():
    beta_report = Path("docs/beta_test_report.md").read_text(encoding="utf-8").lower()
    assumptions = Path("docs/assumptions.md").read_text(encoding="utf-8").lower()

    assert "smoke" in beta_report
    assert "not scientific validation" in beta_report
    assert "not evidence that statistical conclusions are robust" in beta_report
    assert "expected_direction" in assumptions
    assert "multi-gene" in assumptions
    assert "balanced accuracy" in assumptions
    assert "macro-f1" in assumptions


def test_large_sparse_smoke_runs_under_basic_size_budget():
    adata = create_synthetic_perturbseq(n_cells=2500, n_genes=1200, random_state=84)
    adata.X = adata.X.astype("float32")
    split = generate_split(adata, strategy="leave-target-gene-out", random_state=84)

    balance = evaluate_split_balance(adata, split, variables=["batch"])

    assert len(balance) == 1
    assert adata.X.shape == (2500, 1200)
