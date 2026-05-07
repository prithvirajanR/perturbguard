import numpy as np
import anndata as ad
import pandas as pd
from typer.testing import CliRunner

from perturbguard.claims.claim_checker import check_claim
from perturbguard.cli.main import app
from perturbguard.io.config import load_config
from perturbguard.io.schema import apply_schema
from perturbguard.io.split_loader import align_split_to_obs
from perturbguard.leakage.combination_leakage import categorize_combination_leakage
from perturbguard.leakage.split_balance import evaluate_split_balance
from perturbguard.leakage.split_leakage import check_split_leakage
from perturbguard.qc.confounding import detect_confounding, per_perturbation_concentration
from perturbguard.qc.dataset_validator import validate_anndata
from perturbguard.qc.guide_consistency import check_guide_consistency
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.reports.recommendations import build_recommendations
from perturbguard.reports.summary import build_summary
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split


runner = CliRunner()


def test_claim_with_no_test_cells_is_insufficient_metadata():
    adata = create_synthetic_perturbseq(random_state=91)
    split = pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": "train"})

    result = check_claim(adata, split, "unseen_perturbation")

    assert result["status"].iloc[0] == "INSUFFICIENT_METADATA"


def test_align_split_deduplicates_cell_ids_with_last_value_winning():
    adata = create_synthetic_perturbseq(random_state=92)
    cell = str(adata.obs_names[0])
    split = pd.DataFrame(
        {
            "cell_id": [cell, cell],
            "split": ["train", "test"],
        }
    )

    aligned = align_split_to_obs(adata.obs, split)

    assert aligned.loc[adata.obs_names[0]] == "test"


def test_audit_with_missing_perturbation_writes_fail_table_not_traceback(tmp_path):
    adata = create_synthetic_perturbseq(random_state=93)
    adata.obs = adata.obs.drop(columns=["perturbation"])
    data = tmp_path / "missing_perturbation.h5ad"
    out = tmp_path / "audit"
    adata.write_h5ad(data)

    result = runner.invoke(app, ["audit", "--data", str(data), "--out", str(out)])

    assert result.exit_code == 0, result.stdout
    validation = pd.read_csv(out / "tables" / "dataset_validation.csv")
    assert validation.loc[validation["check_id"].eq("obs_perturbation"), "status"].iloc[0] == "FAIL"


def test_guide_consistency_constant_vectors_are_insufficient_not_pass():
    adata = create_synthetic_perturbseq(random_state=94)
    adata.X = np.ones_like(adata.X)

    result = check_guide_consistency(adata)

    assert (result["status"] != "PASS").all()


def test_leakage_checks_preserve_custom_obs_names():
    adata = create_synthetic_perturbseq(random_state=95)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    split = generate_split(adata, strategy="random", random_state=95)

    leakage = check_split_leakage(adata, split, "unseen_perturbation")

    assert not leakage["status"].eq("INSUFFICIENT_METADATA").all()


def test_validation_labels_are_treated_as_evaluation_leakage():
    adata = create_synthetic_perturbseq(random_state=96)
    split = pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": "train"})
    leaked = adata.obs["perturbation"].iloc[0]
    leaked_cells = adata.obs.index[adata.obs["perturbation"].eq(leaked)][:2]
    split.loc[split["cell_id"].isin(leaked_cells.astype(str)), "split"] = "val"

    claim = check_claim(adata, split, "unseen_perturbation")

    assert claim["status"].iloc[0] == "UNSUPPORTED"


def test_all_validation_split_balance_is_insufficient_metadata():
    adata = create_synthetic_perturbseq(random_state=97)
    split = pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": "val"})

    balance = evaluate_split_balance(adata, split)

    assert balance["status"].eq("INSUFFICIENT_METADATA").all()


def test_summary_and_recommendations_treat_unsupported_as_failure():
    sections = {"claim_support": pd.DataFrame([{"status": "UNSUPPORTED"}])}

    assert build_summary(sections)["overall_status"] == "UNSUPPORTED"
    assert not build_recommendations(sections).empty


def test_validate_cli_writes_fail_table_for_corrupted_file(tmp_path):
    bad = tmp_path / "bad.h5ad"
    out = tmp_path / "validation.csv"
    bad.write_bytes(b"not an h5ad")

    result = runner.invoke(app, ["validate", "--data", str(bad), "--out", str(out)])

    assert result.exit_code == 0
    assert pd.read_csv(out)["status"].iloc[0] == "FAIL"


def test_default_config_preserves_boolean_is_control():
    adata = create_synthetic_perturbseq(random_state=98)
    cfg = load_config("configs/default.yaml")

    mapped = apply_schema(adata, cfg.schema, cfg.controls.get("labels"))

    assert mapped.obs["is_control"].sum() == adata.obs["is_control"].sum()


def test_apply_schema_accepts_backed_anndata(tmp_path):
    adata = create_synthetic_perturbseq(random_state=99)
    data = tmp_path / "data.h5ad"
    adata.write_h5ad(data)
    backed = ad.read_h5ad(data, backed="r")

    mapped = apply_schema(backed, {"perturbation": "perturbation"})
    backed.file.close()

    assert mapped.obs.shape[0] == adata.obs.shape[0]


def test_in_memory_split_rejects_invalid_labels():
    adata = create_synthetic_perturbseq(random_state=100)
    split = pd.DataFrame({"cell_id": ["0", "1"], "split": ["train", "holdout"]})

    try:
        align_split_to_obs(adata.obs, split)
    except ValueError as exc:
        assert "Invalid split labels" in str(exc)
    else:
        raise AssertionError("invalid split label was accepted")


def test_empty_anndata_fails_cell_and_gene_count_validation():
    result = validate_anndata(ad.AnnData(np.empty((0, 0))))

    assert result.loc[result["check_id"].eq("n_cells"), "status"].iloc[0] == "FAIL"
    assert result.loc[result["check_id"].eq("n_genes"), "status"].iloc[0] == "FAIL"


def test_all_missing_confounding_metadata_is_insufficient():
    adata = create_synthetic_perturbseq(random_state=101)
    adata.obs["batch"] = pd.NA

    confounding = detect_confounding(adata, metadata_columns=["batch"])
    concentration = per_perturbation_concentration(adata, "batch")

    assert confounding["status"].iloc[0] == "INSUFFICIENT_METADATA"
    assert concentration["status"].eq("INSUFFICIENT_METADATA").all()


def test_signed_expression_target_effect_is_insufficient_metadata():
    adata = create_synthetic_perturbseq(random_state=102)
    adata.X = adata.X - 5

    result = check_target_effect(adata)

    assert result["status"].eq("INSUFFICIENT_METADATA").all()


def test_combination_parser_does_not_split_underscores_by_default_and_splits_hyphen():
    result = categorize_combination_leakage(["Gene_A", "A", "B"], ["Gene_B", "A-B"])

    categories = dict(zip(result["perturbation"], result["category"], strict=False))
    assert categories["Gene_B"] == "STRICT_UNSEEN_COMBINATION"
    assert categories["A-B"] == "RECOMBINATION_OF_SEEN_COMPONENTS"


def test_strict_combination_claim_with_no_combinations_is_insufficient():
    adata = create_synthetic_perturbseq(random_state=103)
    adata = adata[adata.obs["perturbation"].isin(["control", "TP53", "MYC"])].copy()
    split = pd.DataFrame({"cell_id": adata.obs_names.astype(str), "split": "train"})
    split.loc[adata.obs["perturbation"].eq("MYC").values, "split"] = "test"

    result = check_claim(adata, split, "strict_unseen_combination_components")

    assert result["status"].iloc[0] == "INSUFFICIENT_METADATA"


def test_malformed_config_reports_clear_value_error(tmp_path):
    config = tmp_path / "bad.yaml"
    config.write_text("- not\n- a\n- mapping\n", encoding="utf-8")

    try:
        load_config(config)
    except ValueError as exc:
        assert "top-level mapping" in str(exc)
    else:
        raise AssertionError("malformed config was accepted")


def test_underpowered_sparse_confounding_is_not_a_fail():
    adata = create_synthetic_perturbseq(n_cells=120, random_state=104)
    adata.obs["perturbation"] = [f"p{i}" for i in range(adata.n_obs)]

    result = detect_confounding(adata, metadata_columns=["batch"])

    assert result["status"].iloc[0] == "INSUFFICIENT_METADATA"


def test_seen_component_split_holds_out_only_combos_with_seen_components():
    adata = create_synthetic_perturbseq(random_state=105)

    split = generate_split(adata, strategy="seen-component-unseen-combination", random_state=105)
    test_perturbations = set(adata.obs.loc[split["split"].eq("test").values, "perturbation"].astype(str))
    train_perturbations = set(adata.obs.loc[split["split"].eq("train").values, "perturbation"].astype(str))
    train_components = {
        component
        for perturbation in train_perturbations
        for component in perturbation.replace("-", "+").split("+")
        if component and component != "control"
    }

    assert test_perturbations
    assert all("+" in perturbation or "-" in perturbation for perturbation in test_perturbations)
    assert all(set(perturbation.replace("-", "+").split("+")).issubset(train_components) for perturbation in test_perturbations)


def test_seen_component_split_rejects_impossible_dataset():
    adata = create_synthetic_perturbseq(random_state=106)
    adata = adata[adata.obs["perturbation"].isin(["control", "A+B"])].copy()

    try:
        generate_split(adata, strategy="seen-component-unseen-combination", random_state=106)
    except ValueError as exc:
        assert "No seen-component unseen-combination split is possible" in str(exc)
    else:
        raise AssertionError("impossible seen-component split was generated")


def test_split_generation_rejects_duplicate_cell_ids():
    adata = create_synthetic_perturbseq(random_state=107)
    adata.obs_names = ["duplicate"] + [f"cell_{i}" for i in range(adata.n_obs - 1)]
    adata.obs_names = ["duplicate", "duplicate"] + [f"cell_{i}" for i in range(adata.n_obs - 2)]

    try:
        generate_split(adata, strategy="random", random_state=107)
    except ValueError as exc:
        assert "unique observation" in str(exc)
    else:
        raise AssertionError("duplicate cell IDs were accepted")
