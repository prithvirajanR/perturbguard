from pathlib import Path

from perturbguard.reports.recommendations import build_recommendations
from perturbguard.reports.summary import build_summary
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.qc.confounding import detect_confounding
from perturbguard.qc.control_balance import check_control_balance


def test_summary_and_recommendations_are_machine_readable():
    adata = create_synthetic_perturbseq(scenario="batch_confounded", random_state=41)
    sections = {
        "confounding": detect_confounding(adata),
        "control_balance": check_control_balance(adata),
    }

    summary = build_summary(sections)
    recommendations = build_recommendations(sections)

    assert "overall_status" in summary
    assert {"section", "status", "recommendation"}.issubset(recommendations.columns)
    assert not recommendations.empty


def test_norman_demo_script_creates_synthetic_fallback_outputs(tmp_path: Path):
    from scripts.make_demo_subset import main

    out = tmp_path / "norman_demo.h5ad"
    main(["--out", str(out), "--synthetic-fallback"])

    assert out.exists()
