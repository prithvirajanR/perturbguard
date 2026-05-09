# Release Checklist

Before tagging a release:

1. Run `pytest -q`.
2. Run `ruff check .`.
3. Run `python -m build`.
4. Confirm `pyproject.toml`, `perturbguard.__version__`, README release text, and PyPI classifier/maturity language agree.
5. For beta releases, use `Development Status :: 4 - Beta`; do not publish as Production/Stable without the real-data validation below.
6. Run `python scripts/run_smoke_matrix.py --include-real` when real fixtures are present.
7. Review generated reports under `results/smoke_matrix/`.
8. Update `CHANGELOG.md`.
9. Tag the release and publish artifacts.

Production releases should include at least three successful real-data audits:

- SciPlex example
- Jorge example
- LINCS example
