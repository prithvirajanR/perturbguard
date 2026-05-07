from scipy import sparse

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


def test_guide_consistency_is_sparse_safe_and_reports_pairs():
    adata = create_synthetic_perturbseq(scenario="guide_inconsistent", random_state=23)
    adata.X = sparse.csr_matrix(adata.X)

    results = check_guide_consistency(adata)

    assert "n_pairs" in results.columns
    assert (results["n_pairs"] > 0).any()
    assert (results["status"] != "PASS").any()
