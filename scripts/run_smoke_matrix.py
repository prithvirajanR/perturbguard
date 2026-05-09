from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml


def run(cmd: list[str], cwd: Path) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run PerturbGuard CLI smoke matrix.")
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--include-real", action="store_true")
    parser.add_argument("--real-data", default="data/real/sciplex_example.h5ad")
    parser.add_argument("--real-config", default="configs/sciplex_example.yaml")
    args = parser.parse_args(argv)

    cwd = Path(args.workdir).resolve()
    root = cwd / "results" / "smoke_matrix"
    root.mkdir(parents=True, exist_ok=True)

    scenarios = [
        "clean",
        "batch_confounded",
        "low_cell_count",
        "guide_inconsistent",
        "failed_knockdown",
        "leaky_combo",
    ]
    split_strategies = [
        "random",
        "leave-target-gene-out",
        "leave-perturbation-out",
        "strict-unseen-combination",
        "seen-component-unseen-combination",
    ]
    claims = [
        "unseen_perturbation",
        "unseen_target_gene",
        "unseen_combinations",
        "strict_unseen_combination_components",
    ]

    for scenario in scenarios:
        data = root / f"{scenario}.h5ad"
        run(["perturbguard", "simulate", "--scenario", scenario, "--n-cells", "180", "--out", str(data)], cwd)
        run(["perturbguard", "validate", "--data", str(data), "--out", str(root / f"{scenario}_validate.csv")], cwd)
        inferred_config = root / f"{scenario}_inferred.yaml"
        run(["perturbguard", "infer-config", "--data", str(data), "--out", str(inferred_config)], cwd)
        repaired = root / f"{scenario}_repaired.h5ad"
        run(["perturbguard", "repair", "--data", str(data), "--config", str(inferred_config), "--out", str(repaired)], cwd)
        run(["perturbguard", "profile-large", "--data", str(repaired), "--out", str(root / f"{scenario}_profile")], cwd)
        run(["perturbguard", "audit", "--data", str(data), "--out", str(root / f"{scenario}_audit")], cwd)
        run(["perturbguard", "target-map", "--data", str(data), "--out", str(root / f"{scenario}_targets")], cwd)
        run(["perturbguard", "adversarial-check", "--data", str(data), "--out", str(root / f"{scenario}_adversarial")], cwd)
        adata = ad.read_h5ad(data, backed="r")
        predictions = root / f"{scenario}_predictions.csv"
        pd.DataFrame(
            {
                "cell_id": adata.obs_names.astype(str),
                "y_true": adata.obs["perturbation"].astype(str).values,
                "y_pred": adata.obs["perturbation"].astype(str).values,
                "confidence": 0.9,
            }
        ).to_csv(predictions, index=False)
        adata.file.close()
        run(["perturbguard", "evaluate", "--data", str(data), "--predictions", str(predictions), "--out", str(root / f"{scenario}_evaluation")], cwd)
        run(["perturbguard", "dataset-card", "--data", str(data), "--out", str(root / f"{scenario}_dataset_card.md")], cwd)
        for strategy in split_strategies:
            split_out = root / f"{scenario}_{strategy}_split"
            run(
                [
                    "perturbguard",
                    "split",
                    "--data",
                    str(data),
                    "--strategy",
                    strategy,
                    "--out",
                    str(split_out),
                ],
                cwd,
            )
            for claim in claims:
                run(
                    [
                        "perturbguard",
                        "claim",
                        "--data",
                        str(data),
                        "--split",
                        str(split_out / "split.csv"),
                        "--claim",
                        claim,
                        "--out",
                        str(root / f"{scenario}_{strategy}_{claim}_claim"),
                    ],
                    cwd,
                )

    design = root / "design.csv"
    pd.DataFrame(
        {
            "batch": ["b1", "b1", "b2", "b2"],
            "perturbation": ["control", "TP53", "control", "MYC"],
            "n_cells_planned": [100, 100, 100, 100],
        }
    ).to_csv(design, index=False)
    run(["perturbguard", "design-check", "--design", str(design), "--out", str(root / "design")], cwd)
    run(["perturbguard", "power-check", "--design", str(design), "--out", str(root / "power")], cwd)

    run(
        [
            "perturbguard",
            "compare-datasets",
            "--data",
            str(root / "clean.h5ad"),
            "--data",
            str(root / "batch_confounded.h5ad"),
            "--out",
            str(root / "compare"),
        ],
        cwd,
    )
    benchmark_manifest = root / "benchmark.yaml"
    yaml.safe_dump(
        {
            "dataset": str(root / "clean.h5ad"),
            "split": str(root / "clean_leave-perturbation-out_split" / "split.csv"),
            "claim": "unseen_perturbation",
            "model": {"name": "smoke-perfect-baseline"},
            "metrics": ["accuracy", "macro_f1"],
        },
        benchmark_manifest.open("w", encoding="utf-8"),
    )
    run(["perturbguard", "benchmark-check", "--manifest", str(benchmark_manifest), "--out", str(root / "benchmark")], cwd)

    expected_findings = root / "synthetic_expected_findings.csv"
    pd.DataFrame(
        [
            {
                "case": "batch_confounded",
                "section": "metadata_concentration",
                "match_column": "metadata_variable",
                "match_value": "batch",
                "expected_status": "FAIL",
            },
            {
                "case": "failed_knockdown",
                "section": "target_effect",
                "expected_status": "FAIL",
            },
            {
                "case": "low_cell_count",
                "section": "cell_count",
                "expected_status": "FAIL",
            },
        ]
    ).to_csv(expected_findings, index=False)
    validation_manifest = root / "synthetic_validation.yaml"
    yaml.safe_dump(
        {
            "cases": [
                {"name": "batch_confounded", "dataset": str(root / "batch_confounded.h5ad")},
                {"name": "failed_knockdown", "dataset": str(root / "failed_knockdown.h5ad")},
                {"name": "low_cell_count", "dataset": str(root / "low_cell_count.h5ad")},
            ],
            "expected_findings": str(expected_findings),
        },
        validation_manifest.open("w", encoding="utf-8"),
    )
    run(
        [
            "perturbguard",
            "validation-benchmark",
            "--manifest",
            str(validation_manifest),
            "--out",
            str(root / "synthetic_validation"),
            "--fail-on-mismatch",
        ],
        cwd,
    )

    if args.include_real:
        real_cases = [
            {
                "name": "sciplex",
                "data": cwd / args.real_data,
                "config": cwd / args.real_config,
                "split_column": "split_example",
            },
            {
                "name": "jorge",
                "data": cwd / "data/real/Jorge_example.h5ad",
                "config": cwd / "configs/jorge_example.yaml",
                "split_column": "split1",
            },
            {
                "name": "lincs",
                "data": cwd / "data/real/lincs_example.h5ad",
                "config": cwd / "configs/lincs_example.yaml",
                "split_column": "split_example",
            },
        ]
        for case in real_cases:
            if not case["data"].exists():
                print(f"Skipping missing real dataset: {case['data']}", flush=True)
                continue
            run(
                [
                    "perturbguard",
                    "profile-large",
                    "--data",
                    str(case["data"]),
                    "--out",
                    str(root / f"real_{case['name']}_profile"),
                ],
                cwd,
            )
            run(
                [
                    "perturbguard",
                    "audit",
                    "--data",
                    str(case["data"]),
                    "--config",
                    str(case["config"]),
                    "--out",
                    str(root / f"real_{case['name']}_audit"),
                ],
                cwd,
            )
            adata = ad.read_h5ad(case["data"], backed="r")
            real_split = root / f"real_{case['name']}_{case['split_column']}.csv"
            pd.DataFrame(
                {
                    "cell_id": adata.obs_names.astype(str),
                    "split": adata.obs[case["split_column"]].astype(str).values,
                }
            ).to_csv(real_split, index=False)
            adata.file.close()
            run(
                [
                    "perturbguard",
                    "claim",
                    "--data",
                    str(case["data"]),
                    "--split",
                    str(real_split),
                    "--config",
                    str(case["config"]),
                    "--claim",
                    "unseen_perturbation",
                    "--out",
                    str(root / f"real_{case['name']}_claim"),
                ],
                cwd,
            )


if __name__ == "__main__":
    main(sys.argv[1:])
