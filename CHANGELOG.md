# Changelog

## 1.0.0

- First stable PerturbGuard release.
- Adds the full guardrail loop for perturbation benchmarking: config inference, AnnData repair, validation, audit, split generation, claim checking, adversarial leakage checks, model prediction evaluation, benchmark manifest validation, large-file profiling, and dataset-card generation.
- Adds interactive HTML reports with status summary cards, searchable tables, status filters, plot links, and clearer PASS/WARNING/FAIL/SUPPORTED/UNSUPPORTED styling.
- Adds target mapping audits for gene, drug class, pathway/class, missing, and unmapped targets, with optional user-supplied perturbation-to-target maps.
- Adds power/design planning checks for planned cells, controls, replicate support, and batch support.
- Expands CLI smoke testing across synthetic scenarios, real public datasets, split strategies, claim types, design checks, benchmark manifests, and the advanced guardrail commands.
- Packages configs, docs, and Snakemake workflow assets in release artifacts.

## 0.1.0

- Initial production-hardening line for PerturbGuard.
- Adds CLI commands for validation, audit, split generation, claim checking, dataset comparison, design checks, and synthetic simulation.
- Adds schema mapping, synthetic and real-data smoke workflows, split leakage checks, QC tables, recommendations, HTML reports, and plot artifacts.
