# Beta Test Report

PerturbGuard has been smoke-tested on synthetic perturbation scenarios and public example
AnnData files from the PerturbNet reproduction repository.

Real-data cases:

| Dataset | File | Purpose |
|---|---|---|
| SciPlex | `data/real/sciplex_example.h5ad` | Drug perturbation metadata, dose/time/cell-type fields |
| Jorge | `data/real/Jorge_example.h5ad` | Variant perturbation labels and dataset-specific controls |
| LINCS | `data/real/lincs_example.h5ad` | Drug perturbations, many batches/cell IDs, provided holdouts |

The real-data tests verify that schema mapping, missing-control handling, target classification,
audit table writing, HTML reports, plots, and provided split claim checks complete without crashes.
