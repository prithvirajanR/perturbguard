from pathlib import Path

import pandas as pd

from perturbguard.claims.claim_checker import check_claim
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
