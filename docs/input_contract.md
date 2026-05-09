# Input Contract

PerturbGuard accepts AnnData `.h5ad` files and metadata-only design/split tables.

Required canonical observation fields for a full audit:

- `perturbation`
- `is_control`, or perturbation labels from which controls can be inferred

Recommended fields:

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
- `n_counts`
- `pct_mito`

Datasets rarely use these exact names. Use a YAML config to map source columns:

```yaml
schema:
  perturbation: condition
  is_control: is_ctrl
  batch: batch_id
```

Missing optional fields produce `INSUFFICIENT_METADATA` rows rather than crashes.

## Prediction Tables

`perturbguard evaluate` supports two prediction formats.

Class-label prediction rows require:

- `cell_id`
- `y_true`
- `y_pred`

Optional class-label columns:

- `confidence`

Perturbation-effect prediction rows require:

- `perturbation`
- paired `true_<gene>` and `pred_<gene>` columns for observed and predicted mean expression

Optional perturbation-effect columns:

- `pred_low_<gene>` and `pred_high_<gene>` for prediction-interval coverage
- `--gene-sets gene_sets.yaml` where the YAML maps pathway or gene-set names to gene lists

When perturbation-effect columns are present, evaluation writes DE top-k recovery,
effect-vector rank correlation, pathway recovery, and interval calibration tables.

## Validation Benchmark Manifests

Validation benchmark manifests run full audits and compare observed statuses to curated expected
findings.

```yaml
cases:
  - name: batch_confounded
    dataset: data/validation/batch_confounded.h5ad
    config: configs/validation.yaml
    split: splits/batch_confounded.csv
    claim: unseen_perturbation
expected_findings: expected_findings.csv
```

The expected-findings CSV columns are:

- `case`: case name from the manifest
- `section`: audit table name, such as `metadata_concentration` or `target_effect`
- `expected_status`: expected worst status after filtering
- `match_column`: optional column to filter inside the audit table
- `match_value`: optional value to match in `match_column`
