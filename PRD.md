# Molecular Docking Pipeline — Product Requirements Document

**Version:** 1.0
**Author:** Minor Project Team
**Date:** September 2026
**Status:** Foundation (v1.0) — receptor prep, ligand prep, Vina docking, analysis, report

---

## 1. Executive Summary

A self-contained, reproducible **molecular docking pipeline** that takes a protein receptor and a small molecule ligand and produces a ranked set of docking poses with binding-affinity scores and an interpretive report.

The user supplies **either structure files** (PDB receptor + ligand) **or identifiers** (PDB ID / SMILES). The pipeline:

1. Prepares the receptor (PDB → PDBQT).
2. Prepares the ligand (SMILES/SDF → PDBQT).
3. Runs **AutoDock Vina** to sample binding poses.
4. Analyzes results (interaction mapping, affinity ranking, ADMET-lite filter).
5. Emits a **plain-language report** (Markdown + JSON) plus visualization files.

Inspired by the docking module of the open-source [bio-nexus-](https://github.com/Samadsaifi14/bio-nexus-) platform, this project is intentionally scoped **small and self-contained** so it can be run, understood, and extended by a small team or a single student. It is a foundational tool for a drug-discovery / structure-based-design teaching project.

---

## 2. Problem Statement

Molecular docking is a core technique in structure-based drug discovery, yet it is notoriously fragmented to perform:

- Preparing a receptor requires removing water, adding hydrogens, and choosing protonation states.
- Preparing a ligand requires geometry optimization and assignment of Gasteiger charges in a Vina-compatible format.
- Running the docking engine requires a specific search-box configuration.
- Interpreting the raw output requires domain knowledge (what does −8.5 kcal/mol mean? which residues interact?).

A student must juggle PyMOL, ADFR/Meeko, Vina, and analysis scripts, each with its own input/output quirks. This pipeline removes that operational burden: **one command in, one reproducible report out**.

---

## 3. Vision & Goals

**Vision:** A reproducible docking pipeline that any student can run end-to-end without manual file conversion, so they can focus on the biology and interpretation rather than the tooling.

### Goals

| Goal | Metric |
|---|---|
| One-command run | `python -m pipeline run ...` produces an interpretable report |
| Reproducibility | Deterministic with a fixed random seed (Vina `--seed`) |
| Minimal manual steps | No user-side file conversion or PDBQT editing |
| Educational output | Every score is explained in plain language |
| Format coverage | PDB receptor, SMILES/SDF ligand, PDB ID and PubChem ID resolution |

---

## 4. Target Users / Personas

### Primary — Student / Teaching Lab

**Persona: MSc Bioinformatics student**
- Has a protein PDB file (or PDB ID) and a compound (often a known drug, e.g., aspirin, ibuprofen, a kinase inhibitor).
- Knows what docking is but is unfamiliar with the CLI toolchain.
- Needs a correct binding energy and interaction map quickly for a lab report.

### Secondary — Researcher

**Persona: Early-career researcher**
- Wants a quick affinity estimate for a series of compounds against a target prior to expensive screening.
- Needs reproducible results and raw outputs (poses in PDBQT/PDB, affinity table) for downstream use in PyMOL/Chimera.

### Out of Scope (v1)

- Flexible-side-chain docking (rigid receptor only).
- Covalent docking, hydration prediction, water positioning.
- FEP / MM-PBSA binding free energies.
- Explicit-solvent MD refinement.
- Batch virtual-screening of huge libraries (though single multi-ligand batches are supported via CSV).

---

## 5. Pipeline Architecture

```
input: receptor (PDB file | PDB ID)   +   ligand (SMILES | SDF file | compound name)
                │                                │
                ▼                                ▼
        ┌─── RECEPTOR ───┐               ┌─── LIGAND ───┐
        │ sanitize       │               │ parse        │
        │ remove water/  │               │ (RDKit)      │
        │   hetero        │               │ neutralize    │
        │ add H (pH 7.4) │               │ 3D embed      │
        │ assign charges │               │ optimize      │
        │ → receptor.pdbqt│              │ Gasteiger     │
        └───────┬────────┘               └───────┬──────┘
                │                                │
                ▼                                ▼
        ┌───────────── SEARCH BOX ───────────────┐
        │ box center from ligand | pocket        │
        │ (specified coords or ligand centroid)  │
        └──────────────────┬─────────────────────┘
                           ▼
              ┌── AUTODOCK VINA ──┐
              │ stochastic search │
              │ 1.2.5 scoring     │
              │ N modes, seed     │
              └─────────┬─────────┘
                        ▼
      ┌───────────── ANALYSIS ─────────────┐
      │ pose ranking by affinity          │
      │ interaction mapping (H-bond,      │
      │   hydrophobic, π-π, salt bridge)  │
      │ ADMET-lite (Lipinski) from RDKit  │
      └─────────────────┬─────────────────┘
                        ▼
      ┌────────────── REPORT ─────────────┐
      │ results.json (structured)        │
      │ report.md (interpreted)          │
      │ poses.pdbqt / pose_1.pdb ...     │
      │ receptor.pdb + ligand.sdf        │
      └──────────────────────────────────┘
```

### Module responsibilities

| Module | Responsibility | Key tooling |
|---|---|---|
| `receptor_prep` | Normalize receptor PDB → PDBQT (add H, charges) | RDKit / Meeko |
| `ligand_prep` | SMILES/SDF → 3D PDBQT (embed, optimize, charge) | RDKit + Meeko |
| `box` | Compute search-box geometry | geometry / centroid |
| `docking` | Invoke AutoDock Vina subprocess | `bin/vina.exe` |
| `analysis` | Rank poses, map interactions, ADMET-lite | RDKit + custom |
| `report` | Write `results.json`, `report.md`, visuals | stdlib / RDKit |

---

## 6. Input / Output Specification

### Inputs (CLI)

```
python -m pipeline run
    --receptor <pdb_file | PDB_ID>
    --ligand   <smiles | sdf_file | compound_name>
    [--out <dir>]
    [--exhaustiveness 32]      # default 32
    [--num-modes 9]            # default 9
    [--seed 42]                # reproducibility
    [--box-center X,Y,Z]       # optional; default = ligand centroid
    [--box-size 20,20,20]      # Å; optional
```

**Receptor resolution order:** if a file path that exists → use it; else if a 4-char PDB ID → fetch from RCSB; else error.

**Ligand resolution order:** if an SDF/MOL file that exists → use it; else if a comma-free SMILES (parseable by RDKit) → use it; else attempt a PubChem name lookup → else error.

### Outputs (written to `--out`, default `./runs/<timestamp>_<seed>`)

```
results.json         # full structured result (poses, scores, interactions, ADMET)
report.md            # interpreted, human-readable report
poses.pdbqt          # all docking poses in Vina output format
pose_<n>.pdb         # individual pose converted to PDB (for PyMOL/Chimera)
receptor.pdbqt       # prepared receptor
receptor.pdb         # cleaned receptor (H removed-or-added view)
ligand.sdf           # optimized ligand 3D structure
ligand.pdbqt         # prepared ligand
config.txt           # exact Vina config used (reproducibility)
```

---

## 7. Feature Specifications by Phase

### Phase 1 — Core Engine (this PRD baseline) ✅

- **F1.1 Receptor prep:** read PDB/PDB-ID, remove waters & non-protein hetero residues, add hydrogens at pH 7.4, assign Gasteiger charges, emit PDB + PDBQT.
- **F1.2 Ligand prep:** read SMILES/SDF/name, neutralize, embed 3D coordinates, energy-minimize (UFF), assign Gasteiger charges, emit SDF + PDBQT (torsion tree via Meeko).
- **F1.3 Search box:** default box centered on the ligand centroid with configurable size; optional explicit center/size overrides.
- **F1.4 Docking engine:** run AutoDock Vina 1.2.5 with configurable exhaustiveness, num-modes, and fixed seed; capture and parse full log + affinity per mode.
- **F1.5 Pose ranking & output:** rank by best affinity; emit `poses.pdbqt`, per-pose PDB, and structured JSON.
- **F1.6 Interaction mapping:** for the top pose, detect H-bonds, hydrophobic contacts, π–π stacking, and salt bridges against receptor atoms (distance-based, using PDBQT/3D coordinates).
- **F1.7 ADMET-lite:** Lipinski rule-of-five (MW, logP, HBD, HBA) plus top-3 alert flags via RDKit descriptors.
- **F1.8 Interpreted report:** Markdown report that explains each numeric result (affinity bands, best-pose interactions, drug-likeness) in plain language.

### Phase 2 — Depth (planned)

- **F2.1 Pocket-based box:** automated binding-site detection (e.g., fpocket) to auto-center the search box when no cocrystallized ligand exists.
- **F2.2 Multi-ligand screening:** CSV/SDF batch input → ranked table across a ligand set.
- **F2.3 Flexible side chains:** selected receptor side-chain flexibility.
- **F2.4 HTML report:** self-contained interactive HTML (3Dmol.js) alongside Markdown.
- **F2.5 CLI config file:** YAML/JSON config profile instead of many flags.

### Phase 3 — Validation (planned)

- **F3.1 Redocking validation:** redock cocrystallized ligands onto their native receptors; report RMSD against the crystal pose.
- **F3.2 Test suite hardened** to CI; regression-protect interaction detection and parsing.

---

## 8. Non-Functional Requirements

- **Reproducibility:** identical inputs + `--seed` ⇒ byte-comparable pose list and identical rankings.
- **Performance:** a single ligand with default settings completes in ~1–5 minutes on a laptop CPU.
- **Portability:** Python 3.9+ on Windows/Linux/macOS; Vina binary bundled under `bin/`; no GUI required.
- **Discoverability:** outputs are flat, named files in a single run directory; `results.json` is the machine-readable contract (stable schema).
- **Safety:** receptor/ligand are sanitized with RDKit before any geometry/charge work; no remote code; network calls gated to RCSB/PubChem allow-list.
- **Honesty:** the report clearly labels scores as computational predictions, never as measured binding affinities.

---

## 9. Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Receptor prep quality (protonation) | Medium | High | Default neutral-pH handling; document limitations; interactive mode deferred |
| PDB ID unavailable / changed | Medium | Low | Clear error; allow local file fallback; cache downloads |
| PubChem name lookup miss | Medium | Low | Accept SMILES input directly; report unambiguous error |
| Vina binary platform mismatch | Low | High | Bundle per-OS binary under `bin/`; detection + friendly error |
| Interaction mapping false positives | Medium | Medium | Conservative distance thresholds; document method |
| Non-determinism from library versions | Medium | Medium | Pin dependency versions in `requirements.txt` |

---

## 10. Milestone Roadmap

| Milestone | Scope | Status |
|---|---|---|
| M1 — Core engine | Receptor/ligand prep, box, Vina run, parse, JSON | ✅ Shipped |
| M2 — Interpretation | Interaction mapping, ADMET-lite, Markdown report | ✅ Shipped |
| M3 — Validation | Redocking RMSD check, demo datasets, tests | 🔜 In progress |
| M4 — Screening & UI | Multi-ligand batch, HTML/3D report, config profiles | 🔜 Not started |

---

## 11. Open Questions & Future Scope

- **Water & protonation:** should v1 support explicit crystal waters around the binding site? (Deferred — requires per-pdb decision.)
- **Scoring functions:** expose `vinardo` and AD4 as alternatives to the default vina scoring in Phase 2.
- **Validation against experimental data:** build a registry of known complexes to benchmark predicted affinities.
- **Interaction with broader bio-nexus platform:** expose this pipeline as a reusable module/API if the parent platform is later extended.

---

## Appendix A — Environment (as detected on the dev machine)

| Tool | Purpose | Status |
|---|---|---|
| Python 3.14.5 | Runtime | ✅ installed |
| RDKit 2026.3.5 | Cheminformatics (parse/embed/optimize/descriptors) | ✅ installed |
| Meeko 0.8.0 | PDBQT writer (torsion trees) | ✅ installed |
| AutoDock Vina 1.2.5 | Docking engine | ✅ `bin/vina.exe` |
| numpy / scipy / gemmi / biopython | Numerical + structure support | ✅ installed |</parameter>