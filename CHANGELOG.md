# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-09-03

### Added
- **Interactive guided CLI** (`python -m pipeline guided`): a step-by-step
  walk-through that prompts for receptor, ligand, box geometry, and scoring
  settings with validation and sensible defaults.
- **Shared input-validation layer** (`pipeline/validate.py`) reused by both the
  CLI and the web form for consistent, friendly messages.
- **Web form file upload**: receptor PDB, ligand SDF/MOL/PDB, and box-reference
  files can now be uploaded directly in the browser (upload takes precedence
  over typed paths). Invalid uploads/types and conflicting box inputs are
  rejected with clear messages.

### Changed
- Bumped package version to 1.2.0.

## [1.1.1] - 2026-09-03

### Changed
- Receptor PDBQT atom types now written as **element symbols** (`C/N/O/S/P/...`)
  instead of AutoDock4-style names (`OA`/`SA`). Vina receptor scoring is driven by
  element + local chemistry, so this is more correct and reproducible.
- Reports now include a formal **Preparation protocol** section and a clear
  integrity disclaimer (affinities are estimates for ranking, not measurement).

### Security
- The web interface is now **localhost-only by design**: binding to a
  non-loopback host (e.g. `0.0.0.0`) is refused. This matches the intended
  "runs only on my own device" usage and avoids accidental network exposure of
  local file operations.

## [1.1.0] - 2026-09-02

### Added
- FastAPI browser web interface (`python -m pipeline serve`) with a submission
  form, results page, and OpenAPI docs at `/docs`.
- `--box-from` option to center the docking box on a reference pocket file
  (SDF/PDB/PDBQT) in addition to `--box-center` and the ligand centroid default.
- Repository governance documents: `LICENSE` (MIT), `CODE_OF_CONDUCT.md`,
  `CONTRIBUTING.md`, `SECURITY.md` (responsible-use policy), `CHANGELOG.md`,
  and citation metadata (`CITATION.cff`).
- Unit tests for the web interface (`tests/test_web.py`).

### Fixed
- Swapped `rmsd_lb` / `rmsd_ub` column parsing in AutoDock Vina log output.
- Missing `python-multipart` dependency for web form submission.
- Non-ASCII rendering of the dash in the web-server banner on Windows consoles.

## [1.0.0] - 2026-09-01

### Added
- Full 5-stage docking pipeline: ligand prep (Meeko/RDKit), receptor prep,
  search-box definition, AutoDock Vina execution, and analysis + report.
- Outputs: `results.json`, `report.md`, `poses.pdbqt`, per-pose PDB files,
  prepared `receptor.pdbqt` / `ligand.pdbqt`, and `config.txt`.
- Interaction mapping (H-bonds, hydrophobic, salt bridges) and ADMET-lite
  (Lipinski) summaries.
- `pipeline doctor` environment check and reproducibility via `--seed`.
- Demo dataset: streptavidin–biotin (PDB 1STP).