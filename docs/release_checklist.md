# Release Checklist

Before tagging a release:

1. Run `pytest -q`.
2. Run `ruff check .`.
3. Run `python -m build`.
4. Run `python scripts/run_smoke_matrix.py --include-real` when real fixtures are present.
5. Review generated reports under `results/smoke_matrix/`.
6. Update `CHANGELOG.md`.
7. Tag the release and publish artifacts.

Production releases should include at least three successful real-data audits:

- SciPlex example
- Jorge example
- LINCS example
