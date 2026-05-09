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
- Target-effect checking uses configured `expected_direction` when available, infers common
  CRISPRi/knockdown versus CRISPRa/activation directions, and can score configured multi-gene
  targets from `target_genes`. It can flag weak or reversed effects, but it does not prove
  mechanism, specificity, dose response, pathway response, or downstream biological validity.
- Split generation is deterministic and auditable, but it is not a full optimization solver.
  `balanced-random` stratifies over feasible perturbation/metadata groups; rare strata may still
  fall entirely in train when a valid test sample cannot be drawn without creating empty groups.
- Model evaluation supports two modes. Class-label prediction tables report general metrics such
  as balanced accuracy and macro-F1 plus simple confidence calibration. Perturbation-effect
  prediction tables report differential-expression top-k recovery, rank correlation of
  perturbation effects, optional pathway-level recovery from user-provided gene sets, and
  prediction-interval coverage when lower/upper prediction columns are available.
- Target mapping classifies available annotations. Drug mechanisms, pathway annotations, and
  polypharmacology often require curated user-supplied mapping tables.
