from __future__ import annotations

import argparse
from pathlib import Path

from perturbguard.simulate.synthetic_anndata import create_synthetic_perturbseq


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Create a small PerturbGuard demo AnnData file.")
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--synthetic-fallback",
        action="store_true",
        help="Create a synthetic Norman-like fallback instead of requiring a real Norman download.",
    )
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not args.synthetic_fallback:
        raise SystemExit("Real Norman subsetting is not bundled in CI; use --synthetic-fallback.")
    create_synthetic_perturbseq(scenario="clean", n_cells=500, n_genes=80, random_state=2019).write_h5ad(out)


if __name__ == "__main__":
    main()
