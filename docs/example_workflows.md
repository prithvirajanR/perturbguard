# PerturbGuard Example Workflows

## Synthetic Audit

```bash
perturbguard simulate --scenario batch_confounded --out data/demo/batch_confounded.h5ad
perturbguard audit --data data/demo/batch_confounded.h5ad --out results/demo/batch_confounded
```

Open `results/demo/batch_confounded/report.html` to inspect the human-readable report.
Machine-readable outputs are written under `results/demo/batch_confounded/tables/`.

## Generate A Split

```bash
perturbguard split \
  --data data/demo/batch_confounded.h5ad \
  --strategy leave-target-gene-out \
  --out results/demo/split
```

## Check A Claim

```bash
perturbguard claim \
  --data data/demo/batch_confounded.h5ad \
  --split results/demo/split/split.csv \
  --claim unseen_target_gene \
  --out results/demo/claim
```

PerturbGuard reports whether the split supports the requested claim using measurable overlap checks.

## Full Smoke Matrix

```bash
python scripts/run_smoke_matrix.py --include-real
```

This exercises all bundled synthetic scenarios, split strategies, representative claims,
`design-check`, `compare-datasets`, and the downloaded SciPlex example configured by
`configs/sciplex_example.yaml`.

## End-To-End Guardrail Flow

```bash
perturbguard infer-config --data data/raw.h5ad --out configs/inferred.yaml
perturbguard repair --data data/raw.h5ad --config configs/inferred.yaml --out data/repaired.h5ad
perturbguard profile-large --data data/repaired.h5ad --out results/profile
perturbguard audit --data data/repaired.h5ad --config configs/inferred.yaml --out results/audit
perturbguard split --data data/repaired.h5ad --strategy leave-perturbation-out --out results/split
perturbguard claim --data data/repaired.h5ad --split results/split/split.csv --claim unseen_perturbation --out results/claim
perturbguard adversarial-check --data data/repaired.h5ad --split results/split/split.csv --out results/adversarial
perturbguard evaluate --data data/repaired.h5ad --predictions predictions.csv --out results/eval
perturbguard dataset-card --data data/repaired.h5ad --out results/dataset_card.md
```

Prediction files for `evaluate` must include `cell_id`, `y_true`, and `y_pred`.
An optional `confidence` column enables calibration checks.
