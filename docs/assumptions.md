# Statistical Assumptions And Thresholds

PerturbGuard reports measurable audit risks. It does not prove biological truth.

PerturbGuard is guardrail software, not a replacement for study-specific statistical modelling.
The default checks are intended to surface suspicious patterns quickly so that analysts can
inspect them before trusting a benchmark or publication claim.

Default thresholds are intentionally conservative:

- Low perturbation support: warning below 50 cells, fail below 20 cells.
- Batch or metadata confounding: warning at Cramer's V >= 0.3, fail at >= 0.5.
- Dominant metadata concentration: warning at 70%, fail at 90%.
- Metadata-only shortcut baseline: warning at balanced accuracy >= 0.4, fail at >= 0.6.
- Guide consistency: warning at median correlation <= 0.2, fail at <= 0.0.

These are screening thresholds, not universal biological cutoffs. Users should tune them for
experiment size, modality, cell type diversity, dose design, and intended claim.

## Practical Limitations

- Confounding checks use first-pass association tests and effect-size thresholds over available
  metadata. They do not model all nested, paired, longitudinal, dose-response, or mixed-effect
  study designs.
- Split generation is deterministic and auditable, but it is not a full optimization solver.
  `balanced-random` stratifies over feasible perturbation/metadata groups; rare strata may still
  fall entirely in train when a valid test sample cannot be drawn without creating empty groups.
- Model evaluation reports standard classification metrics and simple confidence calibration.
  Serious perturbation prediction benchmarks should add task-specific metrics such as
  differential-expression recovery, rank correlation of perturbation effects, pathway-level
  recovery, and uncertainty/error analysis.
- Target mapping classifies available annotations. Drug mechanisms, pathway annotations, and
  polypharmacology often require curated user-supplied mapping tables.
