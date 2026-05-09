# PerturbGuard

[![CI](https://github.com/prithvirajanR/perturbguard/actions/workflows/ci.yml/badge.svg)](https://github.com/prithvirajanR/perturbguard/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![Release](https://img.shields.io/badge/release-1.0.3-blue)](https://github.com/prithvirajanR/perturbguard/releases/tag/v1.0.3)
[![License](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

**FastQC-style guardrails for single-cell perturbation datasets, splits, claims, and model benchmarks.**

PerturbGuard helps you answer the uncomfortable but necessary question:

> Is this perturbation benchmark actually valid, or did leakage, controls, confounding, weak target effects, or split design make it look better than it is?

It works with AnnData `.h5ad` files and produces interactive HTML reports, machine-readable CSV tables, JSON summaries, and Markdown dataset cards.

**Maturity note:** PerturbGuard is a first public release of guardrail software. Treat it as beta-quality until you have validated its warnings, thresholds, and split behaviour on your own datasets and study designs.

## Why Use It

Perturbation models are easy to overclaim. A model can appear to generalize because the split leaked perturbations, a batch variable predicts the label, controls are missing, target effects failed, or a drug target was treated as a single gene when it is really a pathway/class.

PerturbGuard audits those failure modes before you trust a dataset, split, claim, or benchmark.

## Install

From PyPI:

```bash
pip install perturbguard
```

For local development:

```bash
git clone https://github.com/prithvirajanR/perturbguard.git
cd perturbguard
pip install -e ".[dev]"
pytest -q
```

## Quickstart

```bash
perturbguard simulate --scenario batch_confounded --out data/demo/batch_confounded.h5ad
perturbguard audit --data data/demo/batch_confounded.h5ad --out results/demo_audit
```

Open:

```text
results/demo_audit/report.html
```

You get a searchable report with status filters, summary cards, linked plots, recommendations, and CSV tables under `results/demo_audit/tables/`.

## Common Workflows

### Audit a raw public dataset

```bash
perturbguard infer-config --data data/raw.h5ad --out configs/inferred.yaml
perturbguard repair --data data/raw.h5ad --config configs/inferred.yaml --out data/repaired.h5ad
perturbguard audit --data data/repaired.h5ad --config configs/inferred.yaml --out results/audit
```

### Generate and validate a split

```bash
perturbguard split \
  --data data/repaired.h5ad \
  --strategy leave-perturbation-out \
  --out results/split

perturbguard claim \
  --data data/repaired.h5ad \
  --split results/split/split.csv \
  --claim unseen_perturbation \
  --out results/claim
```

### Evaluate model predictions

```bash
perturbguard evaluate \
  --data data/repaired.h5ad \
  --predictions predictions.csv \
  --out results/evaluation
```

For perturbation-effect predictions, provide perturbation-level expression columns instead of
class labels:

```bash
perturbguard evaluate \
  --data data/repaired.h5ad \
  --predictions effect_predictions.csv \
  --gene-sets gene_sets.yaml \
  --out results/effect_evaluation
```

### Run an expected-findings check

```bash
perturbguard expected-findings \
  --manifest validation.yaml \
  --out results/validation \
  --fail-on-mismatch
```

## Feature Map

| Area | What PerturbGuard Does |
| --- | --- |
| AnnData validation | Checks loadability, matrix presence, cell/gene counts, controls, perturbation metadata, duplicate obs/var names, and malformed files. |
| Config inference | Suggests YAML schema mappings from common metadata aliases in messy public datasets. |
| Repair | Writes normalized `.h5ad` files with canonical metadata, unique indices, inferred controls, and repair action logs. |
| QC audit | Audits perturbation support, controls, cell counts, confounding, metadata shortcuts, target effects, guide consistency, and target mapping. |
| Split generation | Creates random, balanced-random, leave-target-gene-out, leave-perturbation-out, metadata holdout, strict combination, and seen-component combination splits. |
| Claim checking | Verifies whether a split supports claims like unseen perturbation, unseen target gene, or unseen combinations. |
| Leakage detection | Detects train/val/test leakage, combination leakage, unsupported claims, split imbalance, and invalid split labels. |
| Model evaluation | Audits class-label predictions and perturbation-effect predictions, including DE top-k recovery, effect-vector rank correlation, optional pathway recovery, and interval coverage. |
| Adversarial checks | Tests whether metadata-only shortcuts can predict perturbation identity. |
| Target mapping | Classifies targets as measured genes, drug target classes, pathway/class annotations, missing, or unmapped. |
| Design planning | Checks planned cells, controls, replicate support, and batch support before running an experiment. |
| Benchmark manifests | Validates dataset/split/claim/model/metrics manifests and runs claim-support checks. |
| Expected-findings checks | Compares audit outputs against curated expected findings for real or synthetic benchmark cases. |
| Large files | Profiles `.h5ad` files in backed mode before expensive audits. |
| Dataset cards | Generates Markdown dataset cards with audit counts, uses, and limitations. |
| Reports | Writes interactive HTML reports, plot links, CSV tables, `summary.json`, and recommendations. |

## CLI Commands

```text
perturbguard simulate
perturbguard validate
perturbguard audit
perturbguard split
perturbguard claim
perturbguard evaluate
perturbguard repair
perturbguard infer-config
perturbguard target-map
perturbguard compare-datasets
perturbguard design-check
perturbguard power-check
perturbguard benchmark-check
perturbguard expected-findings
perturbguard validation-benchmark
perturbguard profile-large
perturbguard adversarial-check
perturbguard dataset-card
```

## Input Formats

### AnnData

Required for full audits:

- `perturbation`
- `is_control`, or enough configured control labels to infer controls

Recommended:

- `target_gene`
- `guide_id`
- `perturbation_type`
- `batch`
- `replicate`
- `cell_type`
- `donor`
- `plate`
- `timepoint`
- `dose`

### Prediction CSV

For class-label evaluation, required for `perturbguard evaluate`:

- `cell_id`
- `y_true`
- `y_pred`

Optional:

- `confidence`

For perturbation-effect evaluation, use one row per perturbation:

- `perturbation`
- `true_<gene>` columns for observed mean expression
- `pred_<gene>` columns for predicted mean expression
- optional `pred_low_<gene>` and `pred_high_<gene>` columns for interval coverage
- optional `--gene-sets gene_sets.yaml` for pathway-level recovery

### Target Mapping CSV

Recommended for `perturbguard target-map`:

- `perturbation`
- `target`
- optional `target_type`
- optional `source`

### Benchmark Manifest

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

### Expected-Findings Manifest

Expected-findings checks compare curated expected findings against actual audit tables:

```yaml
cases:
  - name: sciplex_known_batch_case
    dataset: data/real/sciplex_case.h5ad
    config: configs/sciplex_example.yaml
expected_findings: validation_expected.csv
```

The expected-findings CSV uses:

- `case`
- `section`
- `expected_status`
- optional `match_column`
- optional `match_value`

## Real Data Smoke Test

PerturbGuard has been smoke-tested locally on public example AnnData files for SciPlex, Jorge, and LINCS-style perturbation datasets.

Example SciPlex run:

```powershell
Invoke-WebRequest `
  -Uri "https://huggingface.co/cyclopeta/PerturbNet_reproduce/resolve/main/example_data/sciplex_example.h5ad" `
  -OutFile data/real/sciplex_example.h5ad

perturbguard audit `
  --data data/real/sciplex_example.h5ad `
  --config configs/sciplex_example.yaml `
  --out results/real/sciplex_audit
```

Full smoke matrix:

```bash
python scripts/run_smoke_matrix.py --include-real
```

## What It Is Not

PerturbGuard is not a perturbation prediction model and does not claim biological truth for unmeasured perturbations. It is a guardrail system: it tells you when your dataset, split, controls, metadata, or benchmark claim may not support the conclusion you want to draw.

Its statistical checks are screening heuristics. Confounding checks, split generation, and model-evaluation metrics are meant to catch common failure modes, not to replace careful study-specific modelling, optimized benchmark design, or perturbation-effect metrics tailored to your biological task.

| PerturbGuard can support | PerturbGuard cannot prove |
| --- | --- |
| Required metadata, controls, split files, and prediction tables are present and parseable. | A perturbation dataset is biologically correct or complete. |
| A supplied split has or does not have observable train/test overlap for configured claims. | A model truly generalizes to all unseen perturbations, targets, doses, or mechanisms. |
| Controls, batches, donors, replicates, and cell counts show measurable warning patterns. | Confounding is absent in unmeasured metadata or hidden experimental structure. |
| Target-effect and model-effect metrics agree with the measured expression fields provided. | A perturbation mechanism, pathway response, or drug polypharmacology is fully validated. |
| Curated expected-findings cases reproduce known audit statuses for tracked datasets. | Warning thresholds are universally calibrated across all real perturbation datasets. |

## Release

Current source version: `1.0.3` beta-quality first public release.

PyPI release metadata should use the Beta classifier from `pyproject.toml` so
the public package metadata matches the maturity language in this repository.

Built and verified with:

- Python 3.11+
- AnnData
- pandas / NumPy / SciPy
- scikit-learn
- Plotly
- Typer
