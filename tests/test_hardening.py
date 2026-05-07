from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from perturbguard.cli.main import app
from perturbguard.io.config import load_config
from perturbguard.io.schema import apply_schema
from perturbguard.io.split_loader import align_split_to_obs
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split
from perturbguard.utils.controls import control_mask


runner = CliRunner()


def test_config_schema_mapping_materializes_canonical_columns(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=61)
    adata.obs["condition_name"] = adata.obs.pop("perturbation")
    adata.obs["plate_id"] = adata.obs.pop("batch")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
schema:
  perturbation: condition_name
  batch: plate_id
  is_control: is_control
""",
        encoding="utf-8",
    )

    mapped = apply_schema(adata, load_config(config_path).schema)

    assert "perturbation" in mapped.obs
    assert "batch" in mapped.obs
    assert mapped.obs["perturbation"].equals(adata.obs["condition_name"])


def test_config_control_labels_materialize_boolean_is_control(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=610)
    adata.obs["condition_name"] = adata.obs["perturbation"]
    adata.obs["control_status"] = adata.obs["perturbation"].where(
        adata.obs["perturbation"].ne("control"),
        "not mutated",
    )
    adata.obs = adata.obs.drop(columns=["is_control", "perturbation"])
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
schema:
  perturbation: condition_name
  is_control: control_status
controls:
  labels: ["not mutated"]
""",
        encoding="utf-8",
    )

    cfg = load_config(config_path)
    mapped = apply_schema(adata, cfg.schema, control_labels=cfg.controls.get("labels"))

    assert mapped.obs["is_control"].dtype == bool
    assert mapped.obs["is_control"].sum() > 0


def test_no_controls_target_effect_returns_insufficient_metadata_not_crash():
    adata = create_synthetic_perturbseq(random_state=62)
    adata.obs["is_control"] = False

    results = check_target_effect(adata)

    assert not results.empty
    assert results["status"].eq("INSUFFICIENT_METADATA").all()


def test_string_false_control_values_are_not_truthy():
    adata = create_synthetic_perturbseq(random_state=620)
    adata.obs["is_control"] = "False"

    assert not control_mask(adata.obs).any()
    results = check_target_effect(adata)
    assert results["status"].eq("INSUFFICIENT_METADATA").all()


def test_split_alignment_uses_cell_ids_and_perturbation_level_splits():
    adata = create_synthetic_perturbseq(random_state=63)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    cell_split = generate_split(adata, strategy="random", random_state=63)

    aligned = align_split_to_obs(adata.obs, cell_split)

    assert len(aligned) == adata.n_obs
    assert set(aligned.unique()) <= {"train", "test", "val"}

    perturbation_split = pd.DataFrame(
        {
            "perturbation": ["control", "TP53", "MYC", "JUN", "FOS", "A+B", "C+D"],
            "split": ["train", "test", "train", "train", "test", "train", "test"],
        }
    )
    aligned = align_split_to_obs(adata.obs, perturbation_split)
    assert aligned.loc[adata.obs["perturbation"].eq("TP53")].eq("test").all()


def test_cli_audit_accepts_config_mapping(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=64)
    adata.obs["condition_name"] = adata.obs.pop("perturbation")
    adata.obs["plate_id"] = adata.obs.pop("batch")
    data = tmp_path / "mapped.h5ad"
    out = tmp_path / "audit"
    config = tmp_path / "config.yaml"
    adata.write_h5ad(data)
    config.write_text(
        """
schema:
  perturbation: condition_name
  is_control: is_control
  target_gene: target_gene
  guide_id: guide_id
  batch: plate_id
  replicate: replicate
  cell_type: cell_type
""",
        encoding="utf-8",
    )

    result = runner.invoke(app, ["audit", "--data", str(data), "--config", str(config), "--out", str(out)])

    assert result.exit_code == 0, result.stdout
    assert (out / "tables" / "dataset_validation.csv").exists()
