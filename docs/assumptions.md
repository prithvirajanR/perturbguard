# Statistical Assumptions And Thresholds

PerturbGuard reports measurable audit risks. It does not prove biological truth.

Default thresholds are intentionally conservative:

- Low perturbation support: warning below 50 cells, fail below 20 cells.
- Batch or metadata confounding: warning at Cramer's V >= 0.3, fail at >= 0.5.
- Dominant metadata concentration: warning at 70%, fail at 90%.
- Metadata-only shortcut baseline: warning at balanced accuracy >= 0.4, fail at >= 0.6.
- Guide consistency: warning at median correlation <= 0.2, fail at <= 0.0.

These are screening thresholds, not universal biological cutoffs. Users should tune them for
experiment size, modality, cell type diversity, dose design, and intended claim.
