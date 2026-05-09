from pathlib import Path

import pandas as pd

from perturbguard.claims.claim_checker import check_claim
from perturbguard.cli.main import _audit_sections
from perturbguard.io.split_loader import load_split
from perturbguard.leakage.combination_leakage import check_combination_leakage
from perturbguard.splitting.strategies import generate_split
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq


def test_split_loader_accepts_cell_level_csv(tmp_path: Path):
    split_path = tmp_path / "split.csv"
    pd.DataFrame({"cell_id": ["0", "1"], "split": ["train", "test"]}).to_csv(split_path, index=False)

    split = load_split(split_path)

    assert split["split"].tolist() == ["train", "test"]
    assert split["cell_id"].tolist() == ["0", "1"]


def test_combination_leakage_table_summarizes_train_test_split():
    adata = create_synthetic_perturbseq(scenario="leaky_combo", random_state=31)
    split = generate_split(adata, strategy="random", test_size=0.3, random_state=31)

    results = check_combination_leakage(adata, split)

    assert {"perturbation", "category", "status"}.issubset(results.columns)
    assert set(results["status"]).issubset({"PASS", "WARNING", "FAIL", "INSUFFICIENT_METADATA"})


def test_claim_checker_supports_donor_batch_and_strict_combination_claims():
    adata = create_synthetic_perturbseq(random_state=32)
    adata.obs["donor"] = adata.obs["batch"].map({"batch1": "donor1", "batch2": "donor2", "batch3": "donor3"})

    donor_split = generate_split(adata, strategy="leave-metadata-out", metadata_column="donor", test_size=0.34)
    donor_claim = check_claim(adata, donor_split, claim="unseen_donor")
    assert donor_claim["status"].iloc[0] == "SUPPORTED"

    combo_split = generate_split(adata, strategy="strict-unseen-combination", test_size=0.3)
    combo_claim = check_claim(adata, combo_split, claim="strict_unseen_combination_components")
    assert combo_claim["status"].iloc[0] in {"SUPPORTED", "INSUFFICIENT_METADATA"}


def test_balanced_random_split_preserves_perturbations_in_train_and_test_when_feasible():
    adata = create_synthetic_perturbseq(n_cells=440, random_state=33)

    split = generate_split(
        adata,
        strategy="balanced-random",
        test_size=0.25,
        random_state=33,
        balance_columns=["batch"],
    )

    merged = adata.obs.assign(split=split["split"].values)
    feasible = merged.groupby(["perturbation", "batch"], observed=True).filter(lambda group: len(group) >= 2)
    counts = feasible.groupby(["perturbation", "batch", "split"], observed=True).size().unstack(fill_value=0)

    assert {"train", "test"}.issubset(set(split["split"]))
    assert (counts.get("test", 0) > 0).any()
    assert (counts.get("train", 0) > 0).all()


def test_audit_with_split_reports_default_claim_support(tmp_path: Path):
    adata = create_synthetic_perturbseq(random_state=34)
    data_path = tmp_path / "data.h5ad"
    split_path = tmp_path / "split.csv"
    adata.write_h5ad(data_path)
    generate_split(adata, strategy="leave-perturbation-out", test_size=0.3).to_csv(split_path, index=False)

    sections = _audit_sections(data_path, split=split_path)

    assert "claim_support" in sections
    assert sections["claim_support"]["claim"].iloc[0] == "unseen_perturbation"
    assert sections["claim_support"]["status"].iloc[0] in {"SUPPORTED", "UNSUPPORTED", "INSUFFICIENT_METADATA"}
