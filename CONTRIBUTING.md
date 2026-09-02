# Contributing to Molecular Docking Pipeline

Thanks for your interest in contributing! This is a small, focused project, so
please read the guidelines below to keep contributions consistent and high
quality.

## Development setup

1. Install Python 3.11+.
2. `pip install -r requirements.txt` (core) and `pip install -r requirements-dev.txt` (tests).
3. Commands:
   - `python -m pipeline doctor` — verify environment
   - `python -m pipeline run ...` — run a docking
   - `python -m pipeline serve` — launch the web UI
   - `python -m pytest tests/ -q` — run the test suite

## Code style

- **No comments unless they clarify non-obvious logic.** Keep code self-documenting.
- Follow PEP 8 with the existing style in `pipeline/` (black-ish, 4-space indent, type hints).
- Draw the line at 100 columns for readability.
- Use `pathlib.Path`, `typing.Optional`/`|`, and the existing logging (`pipeline/utils.py`).

## Workflow

1. **Open an issue first** for anything non-trivial so it can be discussed.
2. Fork the repo, create a branch: `git checkout -b feat/your-feature`.
3. Make focused changes; keep each pull request (PR) small and single-purpose.
4. Add or update tests under `tests/` and run the full suite (must stay green).
5. Update `README.md` if the change affects usage.
6. Submit a PR that closes the linked issue.

## Adding a scientific feature

If you change how structures are prepared, how interactions are computed, or how
scores are interpreted, also update `SECURITY.md` responsible-use notes and the
report "caveats" text so that outputs never overstate predictive accuracy.

## Commit messages

Write commit messages in an imperative, concise style, e.g.:

```
Add clique-based pose clustering (#12)
Fix RMSD column parsing in Vina logs
```

## Licensing

By contributing, you agree that your contributions are licensed under the same
MIT License as the project. See `LICENSE`.