from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData


def create_synthetic_perturbseq(
    scenario: str = "clean",
    n_cells: int = 360,
    n_genes: int = 60,
    random_state: int = 0,
) -> AnnData:
    rng = np.random.default_rng(random_state)
    genes = ["TP53", "MYC", "JUN", "FOS", "A", "B", "C", "D", "E"] + [
        f"Gene{i}" for i in range(max(0, n_genes - 9))
    ]
    genes = genes[:n_genes]
    perturbations = ["control", "TP53", "MYC", "JUN", "FOS", "A", "B", "C", "D", "A+B", "C+D"]
    if scenario == "low_cell_count":
        probs = np.array([0.36, 0.16, 0.14, 0.09, 0.05, 0.04, 0.04, 0.04, 0.04, 0.025, 0.005])
    else:
        probs = np.repeat(1 / len(perturbations), len(perturbations))
    pert = rng.choice(perturbations, size=n_cells, p=probs / probs.sum())
    batches = rng.choice(["batch1", "batch2", "batch3"], size=n_cells)
    plates = rng.choice(["plate1", "plate2", "plate3", "plate4"], size=n_cells)
    if scenario == "batch_confounded":
        by = {
            "control": "batch1",
            "TP53": "batch1",
            "MYC": "batch2",
            "JUN": "batch2",
            "FOS": "batch3",
            "A": "batch1",
            "B": "batch1",
            "C": "batch2",
            "D": "batch2",
            "A+B": "batch3",
            "C+D": "batch3",
        }
        batches = np.array([by[p] for p in pert])
        plates = np.array([f"plate_{p.replace('+', '_')}" for p in pert])
    target = np.array([p.split("+")[0] if p != "control" else "control" for p in pert])
    guide = []
    for p in pert:
        if p == "control":
            guide.append("NTC")
        else:
            guide.append(f"{p.split('+')[0]}_g{rng.integers(1, 3)}")
    obs = pd.DataFrame(
        {
            "perturbation": pert,
            "is_control": pert == "control",
            "target_gene": target,
            "guide_id": guide,
            "perturbation_type": np.where(pert == "control", "control", "CRISPRi"),
            "batch": batches,
            "plate": plates,
            "replicate": rng.choice(["rep1", "rep2"], size=n_cells),
            "cell_type": rng.choice(["K562", "RPE1"], size=n_cells),
        }
    )
    x = rng.poisson(2.0, size=(n_cells, n_genes)).astype(float)
    gene_index = {g: i for i, g in enumerate(genes)}
    for i, p in enumerate(pert):
        if p == "control":
            continue
        for component in p.split("+"):
            if component in gene_index:
                effect = -1.2
                if scenario == "failed_knockdown" and component == "TP53":
                    effect = 1.0
                x[i, gene_index[component]] = np.maximum(0, x[i, gene_index[component]] + effect)
        if scenario == "guide_inconsistent" and obs.iloc[i]["guide_id"].endswith("_g1"):
            x[i, : min(20, n_genes)] = np.maximum(0, x[i, : min(20, n_genes)] - 4)
        if scenario == "guide_inconsistent" and obs.iloc[i]["guide_id"].endswith("_g2"):
            x[i, : min(20, n_genes)] += 5
    adata = AnnData(x, obs=obs, var=pd.DataFrame(index=genes))
    return adata
