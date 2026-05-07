# PerturbGuard

Current release: `1.0.0`

PerturbGuard is a FastQC/MultiQC-style toolkit for single-cell perturbation studies.
It audits perturbation datasets and train/test splits for poor controls, confounding,
low perturbation support, weak guide consistency, failed target effects, and leakage.

PerturbGuard is not a perturbation prediction model and does not claim biological truth
for unmeasured perturbations.

## What it can do

- Validate and repair `.h5ad` AnnData files.
- Infer a starting YAML config from messy public dataset metadata.
- Audit controls, cell-count support, confounding, metadata shortcuts, target effects,
  guide consistency, biological target mapping, and split leakage.
- Generate train/validation/test splits for random, leave-target-gene-out,
  leave-perturbation-out, metadata holdout, and combination-holdout settings.
- Check whether a split supports claims such as unseen perturbation, unseen target gene,
  and unseen combinations.
- Evaluate model prediction files for overall, per-group, and calibration risks.
- Run adversarial leakage checks that ask whether metadata can predict perturbations.
- Profile large `.h5ad` files in backed mode before expensive audits.
- Validate benchmark manifests and planned experimental designs.
- Generate interactive HTML reports and dataset-card Markdown summaries.

## Quickstart

```bash
pip install perturbguard
perturbguard simulate --scenario batch_confounded --out data/demo/batch_confounded.h5ad
perturbguard audit --data data/demo/batch_confounded.h5ad --out results/batch_confounded
```

For local development from this repository:

```bash
pip install -e ".[dev]"
pytest -q
```

## New guardrail commands

```bash
perturbguard infer-config --data data/raw.h5ad --out configs/inferred.yaml
perturbguard repair --data data/raw.h5ad --config configs/inferred.yaml --out data/repaired.h5ad
perturbguard profile-large --data data/repaired.h5ad --out results/profile
perturbguard target-map --data data/repaired.h5ad --mapping targets.csv --out results/targets
perturbguard adversarial-check --data data/repaired.h5ad --out results/adversarial
perturbguard evaluate --data data/repaired.h5ad --predictions predictions.csv --out results/eval
perturbguard power-check --design design.csv --out results/power
perturbguard benchmark-check --manifest benchmark.yaml --out results/benchmark
perturbguard dataset-card --data data/repaired.h5ad --out results/dataset_card.md
```

## Expected input formats

Prediction CSVs for `perturbguard evaluate` must contain:

- `cell_id`
- `y_true`
- `y_pred`

Optional prediction columns:

- `confidence`

Target mapping CSVs for `perturbguard target-map` should contain:

- `perturbation`
- `target`
- optional `target_type`
- optional `source`

Benchmark manifests for `perturbguard benchmark-check` should contain:

```yaml
dataset: data/repaired.h5ad
split: results/split/split.csv
claim: unseen_perturbation
model:
  name: my-model
metrics:
  - accuracy
  - macro_f1
```

## Real-data smoke test

The hardening workflow uses the public SciPlex example AnnData file from Hugging Face:

```bash
Invoke-WebRequest `
  -Uri "https://huggingface.co/cyclopeta/PerturbNet_reproduce/resolve/main/example_data/sciplex_example.h5ad" `
  -OutFile data/real/sciplex_example.h5ad

perturbguard audit `
  --data data/real/sciplex_example.h5ad `
  --config configs/sciplex_example.yaml `
  --out results/real/sciplex_audit
```

To run the CLI smoke matrix across synthetic scenarios, split strategies, claims, design checks,
dataset comparison, and the real SciPlex audit/claim path:

```bash
python scripts/run_smoke_matrix.py --include-real
```
