from perturbguard.qc.cell_count import check_cell_counts
from perturbguard.qc.confounding import per_perturbation_concentration
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq


def test_cell_count_reports_perturbation_batch_and_guide_support():
    adata = create_synthetic_perturbseq(scenario="low_cell_count", random_state=11)

    results = check_cell_counts(
        adata,
        min_cells_warning=50,
        min_cells_fail=20,
        min_cells_per_batch_warning=10,
        min_cells_per_guide_warning=10,
    )

    assert {"entity", "entity_type", "n_cells", "status", "message"}.issubset(results.columns)
    assert {"perturbation", "perturbation_batch", "guide"}.issubset(set(results["entity_type"]))
    assert (results["status"] == "FAIL").any()


def test_per_perturbation_concentration_flags_batch_confounded_data():
    adata = create_synthetic_perturbseq(scenario="batch_confounded", random_state=12)

    results = per_perturbation_concentration(adata, metadata_variable="batch")

    assert {"perturbation", "metadata_variable", "dominant_group", "dominant_fraction", "status"}.issubset(
        results.columns
    )
    assert (results["status"] == "FAIL").any()
