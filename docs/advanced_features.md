# Advanced Features

PerturbGuard now includes ten higher-level guardrails beyond the original audit core.

1. `evaluate`: audits model predictions overall, by metadata group, and by optional confidence calibration.
2. `repair`: writes a normalized AnnData file with canonical metadata, unique indices, and inferred controls where possible.
3. Interactive HTML reports: reports include search, status filters, summary cards, and plot links.
4. `target-map`: classifies target annotations as measured genes, drug classes, pathways/classes, missing, or unmapped.
5. `power-check`: checks planned cells, controls, replicate support, and batch support before an experiment.
6. `benchmark-check`: validates a benchmark manifest and checks whether the declared claim is supported by the split.
7. `profile-large`: opens `.h5ad` files in backed mode and reports shape/storage metadata before expensive audits.
8. `adversarial-check`: asks whether technical metadata can predict perturbation identity.
9. `infer-config`: suggests a YAML schema mapping from common metadata aliases.
10. `dataset-card`: writes a Markdown summary of dataset size, controls, audit statuses, supported uses, and limitations.

These checks are guardrails. They help prevent invalid benchmark claims and brittle analyses, but they do not replace biological validation or model-specific evaluation.
