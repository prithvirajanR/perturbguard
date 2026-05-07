from pathlib import Path

import pandas as pd

from perturbguard.compare.datasets import compare_datasets
from perturbguard.design.checker import check_design
from perturbguard.leakage.leakage_graph import build_leakage_graph
from perturbguard.leakage.pathway_leakage import check_pathway_leakage
from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq
from perturbguard.splitting.strategies import generate_split


def test_design_checker_flags_missing_controls_per_batch():
    design = pd.DataFrame(
        {
            "batch": ["b1", "b1", "b2"],
            "perturbation": ["control", "TP53", "MYC"],
            "n_cells_planned": [100, 100, 100],
        }
    )

    results = check_design(design, min_controls_per_batch=1)

    assert (results["status"] == "FAIL").any()


def test_leakage_graph_and_pathway_leakage_return_structured_tables():
    adata = create_synthetic_perturbseq(random_state=51)
    split = generate_split(adata, strategy="random", test_size=0.25, random_state=51)
    pathways = {"p53_pathway": ["TP53", "MYC"], "ap1_pathway": ["JUN", "FOS"]}

    graph = build_leakage_graph(adata, split)
    pathway = check_pathway_leakage(adata, split, pathways)

    assert {"source", "target", "edge_type"}.issubset(graph.columns)
    assert {"pathway", "train_targets", "test_targets", "status"}.issubset(pathway.columns)


def test_dataset_comparison_and_snakemake_workflow_exist():
    a = create_synthetic_perturbseq(random_state=52)
    b = create_synthetic_perturbseq(random_state=53, n_cells=180)

    comparison = compare_datasets({"a": a, "b": b})

    assert {"dataset", "n_cells", "n_genes", "n_perturbations"}.issubset(comparison.columns)
    assert Path("workflows/snakemake/Snakefile").exists()
