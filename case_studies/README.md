# Case Studies

These case studies are small, versioned examples for checking PerturbGuard against real
AnnData inputs that are already used in the local smoke matrix.

They are expected-findings checks, not scientific proof that the warnings are universally
calibrated. Each case records a dataset path, config path, curated expected findings, and the
expected output of `perturbguard expected-findings` for this release.

Run one case:

```bash
perturbguard expected-findings --manifest case_studies/sciplex/manifest.yaml --out results/case_studies/sciplex --fail-on-mismatch
```

Run all three after downloading the real-data files:

```bash
perturbguard expected-findings --manifest case_studies/sciplex/manifest.yaml --out results/case_studies/sciplex --fail-on-mismatch
perturbguard expected-findings --manifest case_studies/jorge/manifest.yaml --out results/case_studies/jorge --fail-on-mismatch
perturbguard expected-findings --manifest case_studies/lincs/manifest.yaml --out results/case_studies/lincs --fail-on-mismatch
```
