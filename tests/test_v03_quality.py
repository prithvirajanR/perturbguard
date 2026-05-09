from scipy import sparse

from perturbguard.cli.main import _audit_sections
from perturbguard.qc.guide_consistency import check_guide_consistency
from perturbguard.qc.matched_controls import find_matched_controls
from perturbguard.qc.target_effect import check_target_effect
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq


def test_find_matched_controls_prefers_batch_cell_type_timepoint_then_fallback():
    adata = create_synthetic_perturbseq(random_state=21)
    adata.obs["timepoint"] = "24h"
    perturbed_index = adata.obs.index[~adata.obs["is_control"].astype(bool)][0]
    adata.obs.loc[perturbed_index, "batch"] = "unique_batch"

    control_mask, level = find_matched_controls(adata.obs, adata.obs.loc[perturbed_index])

    assert control_mask.any()
    assert level in {"same_cell_type_timepoint", "same_cell_type", "global_controls"}


def test_target_effect_is_sparse_safe_and_reports_matching_level():
    adata = create_synthetic_perturbseq(scenario="failed_knockdown", random_state=22)
    adata.X = sparse.csr_matrix(adata.X)

    results = check_target_effect(adata)

    assert "control_match_level" in results.columns
    assert (results["status"] == "FAIL").any()


def test_target_effect_uses_configured_or_inferred_expected_direction():
    adata = create_synthetic_perturbseq(random_state=221)
    mask = adata.obs["perturbation"].eq("TP53")
    adata.obs.loc[mask, "expected_direction"] = "up"
    target_idx = list(adata.var_names).index("TP53")
    adata.X[mask.values, target_idx] = adata.X[adata.obs["is_control"].astype(bool).values, target_idx].mean() + 4

    configured = check_target_effect(adata)
    tp53 = configured.loc[configured["perturbation"].eq("TP53")].iloc[0]
    assert tp53["expected_direction"] == "up"
    assert tp53["status"] == "PASS"

    adata.obs = adata.obs.drop(columns=["expected_direction"])
    adata.obs.loc[mask, "perturbation_type"] = "CRISPRa"
    inferred = check_target_effect(adata)
    tp53 = inferred.loc[inferred["perturbation"].eq("TP53")].iloc[0]
    assert tp53["expected_direction"] == "up"


def test_target_effect_scores_configured_multi_gene_targets():
    adata = create_synthetic_perturbseq(random_state=222)
    mask = adata.obs["perturbation"].eq("TP53")
    adata.obs.loc[mask, "target_genes"] = "TP53,MYC"
    adata.obs.loc[mask, "expected_direction"] = "any"
    tp53_idx = list(adata.var_names).index("TP53")
    myc_idx = list(adata.var_names).index("MYC")
    adata.X[mask.values, tp53_idx] = adata.X[:, tp53_idx].mean(axis=0) + 5
    adata.X[mask.values, myc_idx] = adata.X[:, myc_idx].mean(axis=0) + 5

    results = check_target_effect(adata)
    row = results.loc[results["perturbation"].eq("TP53")].iloc[0]

    assert row["target_type"] == "multi_gene"
    assert row["n_target_genes"] == 2
    assert row["expected_direction"] == "any"
    assert row["status"] == "PASS"


def test_guide_consistency_is_sparse_safe_and_reports_pairs():
    adata = create_synthetic_perturbseq(scenario="guide_inconsistent", random_state=23)
    adata.X = sparse.csr_matrix(adata.X)

    results = check_guide_consistency(adata)

    assert "n_pairs" in results.columns
    assert (results["n_pairs"] > 0).any()
    assert (results["status"] != "PASS").any()


def test_audit_reports_per_perturbation_metadata_concentration(tmp_path):
    adata = create_synthetic_perturbseq(scenario="batch_confounded", random_state=24)
    path = tmp_path / "confounded.h5ad"
    adata.write_h5ad(path)

    sections = _audit_sections(path)

    assert "metadata_concentration" in sections
    result = sections["metadata_concentration"]
    assert {"perturbation", "metadata_variable", "dominant_fraction", "status"}.issubset(result.columns)
    assert result["metadata_variable"].isin(["batch", "replicate", "donor"]).any()
    assert result["status"].isin(["FAIL", "WARNING"]).any()


def test_audit_preflight_checks_cover_professor_review_points(tmp_path):
    adata = create_synthetic_perturbseq(scenario="batch_confounded", random_state=25)
    path = tmp_path / "audit.h5ad"
    adata.write_h5ad(path)

    sections = _audit_sections(path)

    result = sections["preflight_checks"]
    assert set(result["check_id"]) == {
        "metadata_confounding",
        "split_perturbation_overlap",
        "control_balance",
        "target_effect_detectability",
        "metadata_only_prediction",
        "claim_support",
    }
    assert result.loc[result["check_id"].eq("metadata_confounding"), "status"].iloc[0] == "FAIL"
    assert result.loc[result["check_id"].eq("split_perturbation_overlap"), "status"].iloc[0] == "INSUFFICIENT_METADATA"
    assert result.loc[result["check_id"].eq("claim_support"), "status"].iloc[0] == "INSUFFICIENT_METADATA"
