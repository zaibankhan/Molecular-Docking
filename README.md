# Molecular Docking Pipeline 🧬

A self-contained, reproducible **receptor + ligand → AutoDock Vina docking → interaction report** pipeline. Reference implementation inspired by the docking module of the open-source [bio-nexus](https://github.com/Samadsaifi14/bio-nexus-) platform, scoped as a small standalone project.

One command in → one interpreted report out. No manual PDBQT editing, no GUI, no web tabs.

---

## What it does

```
receptor (PDB file | PDB ID)
        │                           ligand (SMILES | SDF/MOL file | compound name)
        ▼                                        ▼
 prepared receptor.pdbqt            prepared ligand.pdbqt (Meeko torsion tree)
        └──────────────┬─────────────────────────┘
                       ▼
     AutoDock Vina 1.2.5  (search box from ligand centroid, explicit coords,
                           or a reference pocket file)
                       ▼
       pose ranking · interaction map (H-bonds, hydrophobic, salt bridges)
                       ▼
       ADMET-lite (Lipinski) · results.json + report.md + per-pose PDBs
```

The pipeline handles preparation, execution, parsing, and interpretation automatically.

---

## Quick start

```bash
# 1. Install Python deps
pip install -r requirements.txt

# 2. Check the environment (finds vina, rdkit, meeko)
python -m pipeline doctor

# 3. Run a docking
python -m pipeline run \
    --receptor demo/receptor.pdb \
    --ligand "C1C2C(C(=O)CCSC(=O)NCCC1)SCCC2C(=O)O" \
    --box-from demo/pocket_btn.pdb --box-size 15,15,15
```

Output is written to `runs/<timestamp>/`:

```
results.json      structured machine-readable results (stable schema)
report.md         interpreted, plain-language report
poses.pdbqt       all docking poses
poses/pose_N.pdb  individual poses for PyMOL/Chimera
receptor.pbdqt    prepared receptor
ligand.pdbqt      prepared ligand
config.txt        exact Vina invocation (reproducibility)
```

---

## Input options

**Receptor** — a PDB file path *or* a 4-char PDB ID (downloaded from RCSB):

```bash
--receptor 1STP          # fetch from RCSB
--receptor ./receptor.pdb  # local file
```

**Ligand** — auto-resolved in this order: existing `SDF/MOL/PDB` file → valid SMILES → compound name (PubChem lookup):

```bash
--ligand "OC(=O)c1ccccc1C(=O)O"      # phthalic acid by SMILES
--ligand ./demo/ligand.sdf           # by file
--ligand "aspirin"                   # by name (PubChem)
```

**Search box** — three ways to define where Vina samples poses:

| Option | Usage |
|---|---|
| *(default)* | centered on the prepared ligand's centroid |
| `--box-center X,Y,Z` | explicit center (Å) |
| `--box-from FILE` | center on the centroid of an SDF/PDB/PDBQT file (e.g. a cocrystallized ligand or a pocket-defining reference) |
| `--box-size X,Y,Z` | box dimensions in Å (default `20,20,20`) |

---

## Full CLI

```
python -m pipeline doctor                          # env check
python -m pipeline serve [--port 8000]             # local web UI
python -m pipeline run \
    --receptor <path|PDB_ID> \
    --ligand   <smiles|path|name> \
    [--out DIR]            [--exhaustiveness 32]   # search effort
    [--num-modes 9]        [--seed 42]             # poses / reproducibility
    [--ph 7.4]             [-v/--verbose]
    [--box-center X,Y,Z]   [--box-from FILE]       # box placement
    [--box-size 20,20,20]
```

---

## The demo

`demo/` ships a ready-to-run example: the **streptavidin–biotin** complex (PDB `1STP`).

- `demo/receptor.pdb` — streptavidin tetramer receptor
- `demo/pocket_btn.pdb` — the bound biotin coordinates extracted from the crystal, used to center the docking box on the real binding site
- Ligand: biotin by SMILES (or `demo/ligand.sdf`)

```bash
python -m pipeline run \
    --receptor demo/receptor.pdb \
    --ligand "C1C2C(C(=O)CCSC(=O)NCCC1)SCCC2C(=O)O" \
    --box-from demo/pocket_btn.pdb --box-size 15,15,15 \
    --exhaustiveness 24
```

> Note: this is a pedagogical example. Predicted affinities from this simplified rigid-receptor, neutral-pH prep are for ranking and learning, not a substitute for experimental measurements or production-grade protocols.

---

## Web interface (browser-based)

This is a **CLI pipeline** with an optional FastAPI web UI. To run it in your browser on **localhost**:

```bash
python -m pipeline serve --port 8000
```

Then open one of these in your browser:

```
http://127.0.0.1:8000        # main form + results
http://localhost:8000        # same
http://127.0.0.1:8000/docs   # interactive OpenAPI/endpoint docs
```

- `--host 127.0.0.1` (default) binds only to your machine; use `--host 0.0.0.0` to allow other devices on your network.
- `--reload` auto-restarts on code changes (development).
- Web submission writes normal pipeline outputs to `runs/<timestamp>/`; `results.json` and per-pose PDBs are on disk for PyMOL/Chimera.
- Requires the web extras: `pillow`-free, just `fastapi`, `uvicorn[standard]`, `python-multipart`, `httpx` (install via `pip install -r requirements.txt`).

> Note: the web app runs the same `pipeline` code — the browser is just a convenience front-end over the CLI engine.

---

## Running tests

```bash
pip install -r requirements-dev.txt
python -m pytest tests/ -q
```

---

## Reproducibility

Passing `--seed` fixes Vina's random number stream; the same input + seed produces identical poses (verified). `config.txt` records the exact settings for every run.

---

## Project layout

```
├── PRD.md                 # Product requirements document
├── pipeline/              # Python package
│   ├── cli.py             # doctor / run / serve subcommands
│   ├── web/               # FastAPI browser interface (app.py + templates/)
│   ├── orchestrator.py    # wires the stages together
│   ├── receptor_prep.py   # PDB → cleaned PDB + PDBQT
│   ├── ligand_prep.py     # SMILES/SDF → 3D + PDBQT (Meeko/RDKit)
│   ├── docking.py         # Vina subprocess + log/pose parsing
│   ├── analysis.py        # interaction map + ADMET-lite
│   ├── report.py          # results.json + report.md + per-pose PDBs
│   ├── box.py             # search-box geometry
│   ├── deps.py            # environment detection
│   └── utils.py           # PDB/PubChem fetch, run dirs, logging
├── bin/vina.exe           # AutoDock Vina 1.2.5 (Windows)
├── demo/                  # streptavidin–biotin example
├── tests/                 # unit tests
└── requirements*.txt
```

---

## Limitations & honesty

- **Affinities are predictions**, not measurements — useful for *relative* ranking.
- Rigid-receptor docking (no induced fit, no flexible side chains).
- Protonation fixed at the chosen pH; titratable residues may differ in vivo.
- Interaction mapping uses geometric distance thresholds — confirm visually.
- Blind pockets: for high-quality results, center the box on a known site (use `--box-from` or `--box-center`).

## License & attribution

Independent project inspired by the MIT-licensed [`bio-nexus`](https://github.com/Samadsaifi14/bio-nexus-) platform. AutoDock Vina is Apache-2.0 (source: [ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina)) — see the PRD appendix for the environment as bundled.