from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from perturbguard.cli.main import app
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split


runner = CliRunner()


def test_audit_with_split_and_claim_writes_leakage_tables(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=71)
    data = tmp_path / "data.h5ad"
    split_path = tmp_path / "split.csv"
    out = tmp_path / "audit"
    adata.write_h5ad(data)
    generate_split(adata, strategy="leave-target-gene-out", random_state=71).to_csv(split_path, index=False)

    result = runner.invoke(
        app,
        [
            "audit",
            "--data",
            str(data),
            "--split",
            str(split_path),
            "--claim",
            "unseen_target_gene",
            "--out",
            str(out),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert (out / "tables" / "split_leakage.csv").exists()
    assert (out / "tables" / "combination_leakage.csv").exists()
    assert (out / "tables" / "claim_support.csv").exists()


def test_compare_datasets_and_design_check_cli(tmp_path: Path):
    a = tmp_path / "a.h5ad"
    b = tmp_path / "b.h5ad"
    compare_out = tmp_path / "compare"
    design_path = tmp_path / "design.csv"
    design_out = tmp_path / "design_report"
    create_synthetic_perturbseq(random_state=72).write_h5ad(a)
    create_synthetic_perturbseq(random_state=73, n_cells=120).write_h5ad(b)
    pd.DataFrame(
        {
            "batch": ["b1", "b1", "b2"],
            "perturbation": ["control", "TP53", "MYC"],
            "n_cells_planned": [100, 100, 20],
        }
    ).to_csv(design_path, index=False)

    result = runner.invoke(app, ["compare-datasets", "--data", str(a), "--data", str(b), "--out", str(compare_out)])
    assert result.exit_code == 0, result.stdout
    assert (compare_out / "dataset_comparison.csv").exists()

    result = runner.invoke(app, ["design-check", "--design", str(design_path), "--out", str(design_out)])
    assert result.exit_code == 0, result.stdout
    assert (design_out / "design_validation.csv").exists()
    assert (design_out / "report.html").exists()


def test_claim_cli_accepts_config_mapping(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=74)
    adata.obs["condition_name"] = adata.obs.pop("perturbation")
    data = tmp_path / "mapped.h5ad"
    split = tmp_path / "split.csv"
    config = tmp_path / "config.yaml"
    out = tmp_path / "claim"
    adata.write_h5ad(data)
    pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": ["train"] * 300 + ["test"] * 60}).to_csv(
        split,
        index=False,
    )
    config.write_text(
        """
schema:
  perturbation: condition_name
  is_control: is_control
  target_gene: target_gene
""",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "claim",
            "--data",
            str(data),
            "--split",
            str(split),
            "--config",
            str(config),
            "--claim",
            "unseen_perturbation",
            "--out",
            str(out),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert (out / "claim_support.csv").exists()


def test_audit_report_links_generated_plots_and_fail_on_option(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=75)
    adata.obs = adata.obs.drop(columns=["perturbation"])
    data = tmp_path / "bad.h5ad"
    out = tmp_path / "audit"
    adata.write_h5ad(data)

    result = runner.invoke(
        app,
        ["audit", "--data", str(data), "--out", str(out), "--fail-on", "fail"],
    )

    assert result.exit_code == 1
    assert (out / "summary.json").exists()

    good = tmp_path / "good.h5ad"
    good_out = tmp_path / "good_audit"
    create_synthetic_perturbseq(random_state=76).write_h5ad(good)
    result = runner.invoke(app, ["audit", "--data", str(good), "--out", str(good_out)])
    assert result.exit_code == 0, result.stdout
    html = (good_out / "report.html").read_text(encoding="utf-8")
    assert "plots/" in html


def test_split_cli_exposes_seed_and_validation_size(tmp_path: Path):
    data = tmp_path / "data.h5ad"
    out1 = tmp_path / "split1"
    out2 = tmp_path / "split2"
    create_synthetic_perturbseq(random_state=77).write_h5ad(data)

    result = runner.invoke(
        app,
        [
            "split",
            "--data",
            str(data),
            "--out",
            str(out1),
            "--strategy",
            "random",
            "--random-state",
            "1",
            "--val-size",
            "0.1",
        ],
    )
    assert result.exit_code == 0, result.stdout
    result = runner.invoke(
        app,
        [
            "split",
            "--data",
            str(data),
            "--out",
            str(out2),
            "--strategy",
            "random",
            "--random-state",
            "2",
        ],
    )
    assert result.exit_code == 0, result.stdout

    split1 = pd.read_csv(out1 / "split.csv")
    split2 = pd.read_csv(out2 / "split.csv")
    assert "val" in set(split1["split"])
    assert not split1["split"].equals(split2["split"])
