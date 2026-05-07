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
