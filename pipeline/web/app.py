"""FastAPI web interface for the molecular docking pipeline.

Purpose: a clean, reliable browser front-end to the CLI/pipeline engine.
The user fills in receptor + ligand + box, clicks "Run docking", and is shown
a progress page while Vina works, then lands on a results page.

Run with:  python -m pipeline serve
Then open http://127.0.0.1:8000 in your browser. Localhost-only by design.
"""
from __future__ import annotations

import datetime as _dt
import json
import threading
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from ..config import DockingConfig
from .. import validate as vld
from .. import seqlab
from ..utils import fetch_pdb, fetch_rcsb_entry, looks_like_pdb_id, LOG

_HERE = Path(__file__).resolve().parent
_TEMPLATES_DIR = _HERE / "templates"
_STATIC_DIR = _HERE / "static"
_RUNS_BASE = Path("runs")

dapp = FastAPI(title="Molecular Docking Pipeline", version="1.4.0")
dapp.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")
templates = Jinja2Templates(directory=str(_TEMPLATES_DIR))


def _aa_color(ch: str) -> str:
    """Clustal-style amino-acid colouring used by the MSA result view."""
    if not ch or ch == "-":
        return "transparent"
    p = ch.upper()
    if p in "ILVMAFW":
        return "#0057ae"
    if p in "GSTYC":
        return "#00b000"
    if p in "DENQ":
        return "#c00000"
    if p in "KRH":
        return "#6f00d0"
    if p == "P":
        return "#8f4700"
    return "#000000"


def _ws_min(a, b):
    return a if a < b else b


templates.env.globals["color"] = _aa_color
templates.env.globals["ws_min"] = _ws_min

# Shared in-memory state for the most recent run (thread-safe).
_STATE: dict = {"running": False, "results": None, "out_dir": None, "error": None}
_LOCK = threading.Lock()

_MAX_UPLOAD = 100 * 1024 * 1024  # 100 MB


def _render(request: Request, template: str, **ctx):
    return templates.TemplateResponse(request, template, ctx)


@dapp.get("/", response_class=HTMLResponse)
def index(request: Request):
    with _LOCK:
        running = _STATE["running"]
    return _render(request, "index.html", running=running)


@dapp.get("/dock", response_class=HTMLResponse)
def dock_page(request: Request, running: bool = False):
    with _LOCK:
        _running = _STATE["running"]
    return _render(request, "dock.html", running=running or _running)


@dapp.get("/api/health")
def health():
    return {"status": "ok"}


@dapp.get("/api/status")
def status():
    with _LOCK:
        return {
            "running": _STATE["running"],
            "has_results": _STATE["results"] is not None,
            "error": _STATE["error"],
            "out_dir": _STATE["out_dir"],
        }


@dapp.get("/api/resolve_ligand")
def api_resolve_ligand(text: str = ""):
    """Resolve a ligand input via PubChem: file/SMILES pass through; a compound
    name (e.g. 'aspirin') is converted to its canonical SMILES with properties."""
    try:
        resolved, name, record = _resolve_ligand_input(text)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    result = {
        "ok": True,
        "smiles": resolved,
        "name": name,
        "action": "resolved" if name else "direct",
    }
    if record:
        result["properties"] = {
            "cid": record.get("cid"),
            "formula": record.get("formula"),
            "molecular_weight": record.get("molecular_weight"),
            "iupac_name": record.get("iupac_name"),
        }
    return result


@dapp.get("/api/resolve_receptor")
def api_resolve_receptor(pdb: str = ""):
    """Fetch RCSB metadata for a 4-character PDB ID (title, resolution, method)."""
    pdb_id = (pdb or "").strip().upper()
    if not looks_like_pdb_id(pdb_id):
        return {"ok": False, "error": "Enter a valid 4-character PDB ID, e.g. 1STP."}
    try:
        info = fetch_rcsb_entry(pdb_id)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": str(exc)}
    return {"ok": True, **info}


def _run_job(cfg: DockingConfig) -> None:
    """Run the pipeline in a background thread and publish the outcome to _STATE."""
    from ..orchestrator import run

    with _LOCK:
        _STATE["running"] = True
        _STATE["error"] = None
        _STATE["results"] = None
        _STATE["out_dir"] = None
    try:
        results = run(cfg, verbose=False)
        out_dir = cfg.resolved_out()
        with _LOCK:
            _STATE["results"] = results
            _STATE["out_dir"] = str(out_dir)
    except Exception as exc:  # noqa: BLE001
        with _LOCK:
            _STATE["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        with _LOCK:
            _STATE["running"] = False


async def _persist_upload(upload: Optional[UploadFile], name: str, allowed: tuple,
                          staging: Path) -> Optional[str]:
    """Write an uploaded file to the staging dir. Returns its path or None.

    Extension-less uploads are sniffed by content (PDB / PDBQT / SDF markers)
    and saved with the correct extension. If the content cannot be classified,
    the file is still accepted with a sensible default extension for the slot,
    so the pipeline itself can report the real format problem downstream.
    Only a mismatched *declared* extension is rejected outright.
    """
    if upload is None:
        return None
    filename = upload.filename or ""
    ext = Path(filename).suffix.lower()
    data = await upload.read()
    if not data:
        LOG.warning("Uploaded %s file is empty; ignoring it.", name)
        return None
    if len(data) > _MAX_UPLOAD:
        raise ValueError(f"Uploaded {name} file is too large (>100 MB).")
    if ext not in allowed:
        probe = data[:4096].decode("utf-8", errors="replace")
        sniffed = None
        if any(m in probe for m in ("REMARK VINA RESULT", "ROOT", "BRANCH", "TORSDOF", "@<TRIPOS>")):
            sniffed = ".pdbqt"
        elif any(m in probe for m in ("HEADER", "CRYST1", "CONECT", "ATOM  ", "HETATM", "TER   ")):
            sniffed = ".pdb"
        elif any(m in probe for m in ("M  END", "$$$$", "V2000", "V3000")):
            sniffed = ".sdf"
        elif ".smi" in allowed:
            first = next((ln.strip() for ln in probe.splitlines() if ln.strip()), "")
            if (first and len(first) < 250
                    and all(c not in " \t" for c in first) and any(c.isalpha() for c in first)):
                sniffed = ".smi"
        if sniffed and sniffed in allowed:
            ext = sniffed
        elif not ext:
            default_ext = {"receptor": ".pdb", "ligand": ".sdf", "box": ".pdb"}.get(name, ".pdb")
            ext = default_ext if default_ext in allowed else allowed[0]
        else:
            raise ValueError(
                f"Unsupported file type for {name}: '{ext or 'no-extension'}'. "
                f"Allowed: {', '.join(allowed)}"
            )
    dest = staging / f"{name}_{staging_token()}{ext}"
    dest.write_bytes(data)
    return str(dest)


def _looks_like_structure(text: str) -> bool:
    """True if text contains markers of a molecular structure file (SDF/MOL,
    PDB, PDBQT/MOL2) rather than a bare SMILES line."""
    markers = (
        "HEADER", "CRYST1", "CONECT", "ATOM  ", "HETATM", "TER   ",
        "M  END", "$$$$", "V2000", "V3000",
        "REMARK VINA RESULT", "ROOT", "BRANCH", "TORSDOF", "@<TRIPOS>",
    )
    return any(m in text for m in markers)


def _read_smiles_from_file(path: Path) -> Optional[str]:
    """Return the first non-empty token of a SMILES text file, or None."""
    try:
        payload = path.read_text(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        return None
    if _looks_like_structure(payload[:4096]):
        return None
    for line in payload.splitlines():
        line = line.strip()
        if line:
            return line.split()[0] if line.split()[0] else line
    return None


def staging_token() -> str:
    import time

    return str(int(time.time() * 1000))


def _resolve_ligand_input(value: str) -> tuple[str, Optional[str], Optional[dict]]:
    """Return (ligand_value, source_name, pubchem_record).

    ``value`` may be an existing file path or a valid SMILES (both kept as-is,
    ``source_name=None``, ``pubchem_record=None``) or a compound name that is
    resolved via the PubChem REST API (``source_name`` = the resolved name and
    ``pubchem_record`` = its SMILES/formula/properties). Raises ValueError if
    none of the three interpretations work.
    """
    value = (value or "").strip()
    if not value:
        return value, None, None
    if Path(value).exists():
        return value, None, None
    from rdkit import Chem

    if Chem.MolFromSmiles(value):
        return value, None, None
    from ..utils import fetch_pubchem_record

    try:
        record = fetch_pubchem_record(value)
    except Exception as exc:  # noqa: BLE001
        raise ValueError(
            f"Could not resolve ligand '{value}': it is not a file, not a valid SMILES, "
            f"and the PubChem lookup failed ({type(exc).__name__}: {exc})."
        ) from exc
    return record["smiles"], value, record


async def _start_docking(
    request: Request,
    *,
    receptor: str,
    ligand: str,
    box_center: str = "",
    box_source: str = "",
    box_size: str = "20,20,20",
    exhaustiveness: int = 16,
    num_modes: int = 9,
    seed: int = 42,
    ph: float = 7.4,
    receptor_file: Optional[UploadFile] = None,
    ligand_file: Optional[UploadFile] = None,
    ligand_smiles_file: Optional[UploadFile] = None,
    box_file: Optional[UploadFile] = None,
    error_template: str = "dock.html",
):
    """Validate inputs, persist uploads, and launch a background docking run."""
    errors: list[str] = []

    receptor = (receptor or "").strip()
    ligand = (ligand or "").strip()
    box_center = (box_center or "").strip()
    box_source = (box_source or "").strip()
    box_size = (box_size or "").strip()

    # An uploaded SMILES text file takes priority over the ligand text field.
    if ligand_smiles_file is not None:
        try:
            payload = (await ligand_smiles_file.read()).decode("utf-8", errors="replace")
            for line in payload.splitlines():
                line = line.strip()
                if line:
                    ligand = line.split()[0] if line.split()[0] else line
                    break
        except Exception:  # noqa: BLE001
            ligand = ligand

    # Resolve a compound name (e.g. "aspirin") to its canonical SMILES via PubChem
    # before launching, so a bad name fails fast instead of failing mid-run.
    ligand_as_typed = ligand
    resolved_note = None
    if ligand:
        try:
            ligand, resolved_name, _ = _resolve_ligand_input(ligand)
        except ValueError as exc:
            errors.append(str(exc))
        else:
            if resolved_name:
                resolved_note = (
                    f"Resolved '{resolved_name}' via PubChem — docking the canonical "
                    f"SMILES: {ligand}"
                )

    # Download a receptor PDB ID from RCSB eagerly so a missing ID fails fast
    # while the user is on the form, and the background job needs no network.
    pdb_id_path = None
    if receptor and receptor_file is None and looks_like_pdb_id(receptor):
        try:
            _RUNS_BASE.mkdir(parents=True, exist_ok=True)
            staging = _RUNS_BASE / "_uploads"
            staging.mkdir(parents=True, exist_ok=True)
            dest = staging / f"receptor_{receptor.lower()}.pdb"
            pdb_id_path = str(fetch_pdb(receptor.upper(), dest))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"Could not download receptor PDB '{receptor}': {exc}")

    if not receptor and receptor_file is None:
        errors.append("Receptor is required (upload a PDB file or type a path / PDB ID).")
    if not ligand and ligand_file is None:
        errors.append(
            "Ligand is required (upload an SDF/MOL/PDB, a SMILES (.smi/.txt) file, "
            "or type SMILES / path / name)."
        )

    try:
        box_center_tuple = vld.parse_triple(box_center, "box center")
    except vld.InputError as exc:
        errors.append(str(exc))
        box_center_tuple = None
    try:
        box_size_tuple = vld.parse_triple(box_size, "box size") or (20.0, 20.0, 20.0)
    except vld.InputError as exc:
        errors.append(str(exc))
        box_size_tuple = (20.0, 20.0, 20.0)

    if box_center and box_source:
        errors.append("Provide EITHER a box center OR a box-source file, not both.")

    prev = {"receptor": receptor, "ligand": ligand_as_typed, "box_center": box_center,
            "box_source": box_source, "box_size": box_size,
            "exhaustiveness": exhaustiveness, "num_modes": num_modes, "seed": seed, "ph": ph}

    if errors:
        return _render(request, error_template, running=False, errors=errors, prev=prev)

    # --- Persist uploaded files. ---
    staging = _RUNS_BASE / "_uploads"
    staging.mkdir(parents=True, exist_ok=True)
    try:
        rec_path = await _persist_upload(receptor_file, "receptor", (".pdb", ".ent", ".pdbqt"), staging)
        lig_path = await _persist_upload(
            ligand_file, "ligand", (".sdf", ".mol", ".pdb", ".pdbqt", ".smi", ".txt"), staging
        )
        box_path = await _persist_upload(box_file, "box", (".sdf", ".pdb", ".pdbqt"), staging)
    except ValueError as exc:
        return _render(request, error_template, running=False, errors=[str(exc)], prev=prev)

    # A ligand file whose content is not a structure is a SMILES text file:
    # use its first token as the ligand (works for .smi/.txt and for
    # extension-less files that were saved with a default structure extension).
    smi_note = None
    if lig_path:
        p = Path(lig_path)
        smi = _read_smiles_from_file(p)
        if smi is not None:
            ligand = smi
            lig_path = None
            smi_note = f"Using SMILES '{smi}' from the uploaded file."

    receptor_value = rec_path or pdb_id_path or receptor
    ligand_value = lig_path or ligand

    if not ligand_value.strip():
        return _render(
            request, error_template, running=False,
            errors=["Ligand is required (upload a non-empty SDF/MOL/PDB, a SMILES "
                    "(.smi/.txt) file, or type SMILES / path / name)."],
            prev=prev,
        )
    if not receptor_value.strip():
        return _render(
            request, error_template, running=False,
            errors=["Receptor is required (upload a non-empty PDB file or type a "
                    "PDB ID / path)."],
            prev=prev,
        )
    box_source_value = box_path or (box_source or None)

    cfg = DockingConfig(
        receptor=receptor_value,
        ligand=ligand_value,
        out_dir=None,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        ph=ph,
        box_center=box_center_tuple,
        box_size=box_size_tuple,
        box_source=box_source_value,
    )

    # Launch in the background; the browser is redirected to the running page.
    with _LOCK:
        _STATE["running"] = True
    threading.Thread(target=_run_job, args=(cfg,), daemon=True).start()
    note_parts = [n for n in (resolved_note, smi_note) if n]
    return _render(request, "running.html", running=True, note=" · ".join(note_parts) or None)


@dapp.post("/submit")
async def submit(
    request: Request,
    receptor: str = Form(""),
    ligand: str = Form(""),
    box_center: str = Form(""),
    box_source: str = Form(""),
    box_size: str = Form("20,20,20"),
    exhaustiveness: int = Form(16),
    num_modes: int = Form(9),
    seed: int = Form(42),
    ph: float = Form(7.4),
    receptor_file: Optional[UploadFile] = File(None),
    ligand_file: Optional[UploadFile] = File(None),
    ligand_smiles_file: Optional[UploadFile] = File(None),
    box_file: Optional[UploadFile] = File(None),
):
    return await _start_docking(
        request,
        receptor=receptor,
        ligand=ligand,
        box_center=box_center,
        box_source=box_source,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        ph=ph,
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        ligand_smiles_file=ligand_smiles_file,
        box_file=box_file,
    )


@dapp.get("/results", response_class=HTMLResponse)
def results(request: Request, out: Optional[str] = None):
    with _LOCK:
        error = _STATE["error"]
        in_mem = _STATE["results"]
        in_mem_out = _STATE["out_dir"]
        running = _STATE["running"]

    out_dir = out or in_mem_out
    results_dict = None
    if out_dir:
        results_json = Path(out_dir) / "results.json"
        if results_json.exists():
            try:
                results_dict = json.loads(results_json.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                results_dict = None
    if results_dict is None:
        results_dict = in_mem

    if running and results_dict is None and error is None:
        return _render(request, "running.html", running=True)
    if error and results_dict is None:
        return _render(request, "results.html", results=None, out_dir=out_dir, error=error)
    return _render(request, "results.html", results=results_dict, out_dir=out_dir, error=None)


# ---------------------------------------------------------------------------
# History: catalog past runs saved under the runs/ directory.
# ---------------------------------------------------------------------------

def _parse_run_ts(name: str) -> str:
    """Parse a YYYYMMDD_HHMMSS run-dir name into a friendly local timestamp."""
    try:
        dt = _dt.datetime.strptime(name, "%Y%m%d_%H%M%S")
        return dt.strftime("%d %b %Y, %H:%M:%S")
    except ValueError:
        return name


def _list_runs(limit: Optional[int] = None) -> list[dict]:
    """Scan runs/ and return metadata for every completed run, newest first."""
    runs_dir = _RUNS_BASE.resolve()
    entries = []
    if not runs_dir.is_dir():
        return entries
    for folder in runs_dir.iterdir():
        if not folder.is_dir() or folder.name.startswith("_"):
            continue
        rj = folder / "results.json"
        if not rj.is_file():
            continue
        try:
            data = json.loads(rj.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            continue
        docking = data.get("docking") or {}
        best = docking.get("best_pose") or {}
        # Skip incomplete/error runs that only have a bare docking stub.
        if "schema_version" not in data:
            continue
        ligand = data.get("ligand") or {}
        receptor = data.get("receptor") or {}
        cfg = data.get("config") or {}
        entries.append({
            "name": folder.name,
            "out_dir": str(folder),
            "run_ts": _parse_run_ts(folder.name),
            "receptor": receptor.get("source", "—"),
            "ligand_input": ligand.get("input", "—"),
            "formula": ligand.get("formula", "—"),
            "mw": ligand.get("mol_weight", "—"),
            "affinity": best.get("affinity_kcal_mol", "—"),
            "num_poses": docking.get("num_modes", docking.get("pose_files", 0)),
            "pose_files": len(data.get("pose_files") or []),
            "box_size": cfg.get("box_size", "—"),
            "report_exists": (folder / "report.md").is_file(),
        })
    entries.sort(key=lambda e: e["name"], reverse=True)
    if limit:
        entries = entries[:limit]
    return entries


_SEQLAB_HISTORY = _RUNS_BASE / "_seqlab_history.json"
_SEQLAB_HISTORY_LIMIT = 40


def _load_seqlab_history() -> list[dict]:
    p = _SEQLAB_HISTORY.resolve()
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return []
    return data if isinstance(data, list) else []


def _save_seqlab_history(records: list[dict]) -> None:
    try:
        _SEQLAB_HISTORY.resolve().parent.mkdir(parents=True, exist_ok=True)
        _SEQLAB_HISTORY.write_text(
            json.dumps(records, ensure_ascii=False, indent=1), encoding="utf-8")
    except Exception as exc:  # noqa: BLE001
        LOG.warning("could not save sequence-lab history: %s", exc)


def _seq_key() -> str:
    """Unique, sortable key for a sequence-lab record (microsecond resolution)."""
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")


def _append_seqlab(record: dict) -> None:
    records = _load_seqlab_history()
    records.insert(0, record)
    _save_seqlab_history(records[:_SEQLAB_HISTORY_LIMIT])


def _seqlab_rows() -> list[dict]:
    """Build display rows for the history page from saved BLAST/MSA records."""
    rows = []
    for r in _load_seqlab_history():
        ts = r.get("ts", "")
        display = _parse_run_ts(ts[:15]) if len(ts) >= 15 else ts
        if r.get("kind") == "blast":
            hits = r.get("hits") or []
            q = (r.get("query") or "").strip()[:40]
            rows.append({
                "detail": ts,
                "run_ts": display,
                "kind": "blast",
                "label": "BLAST",
                "summary": (f"query '{q}' -> {len(hits)} hit(s) | "
                            f"{r.get('db_size', '?')} subject(s) in database"),
            })
        else:
            seqs = r.get("sequences") or []
            rows.append({
                "detail": ts,
                "run_ts": display,
                "kind": "msa",
                "label": "MSA",
                "summary": (f"{len(seqs)} sequences | {r.get('columns', '?')} columns | "
                            f"{r.get('method', '')}"),
            })
    return rows


@dapp.get("/history", response_class=HTMLResponse)
def history_page(request: Request):
    return _render(request, "history.html", runs=_list_runs(),
                   seqlab=_seqlab_rows())


@dapp.get("/seqlab", response_class=HTMLResponse)
def seqlab_view(request: Request, ts: str = ""):
    """Re-render a saved BLAST search or MSA from sequence-lab history."""
    for rec in _load_seqlab_history():
        if rec.get("ts") != ts:
            continue
        if rec.get("kind") == "blast":
            hits = rec.get("hits") or []
            lines = [f">query\n{rec.get('query', '')}"]
            for h in hits[:10]:
                seq = h.get("subject_seq") or ""
                if seq:
                    lines.append(f">{h['subject_id']}\n{seq}")
            return _render(request, "blast.html", results=hits, error=None,
                           query=rec.get("query", ""), db_size=rec.get("db_size", 0),
                           sensitive=bool(rec.get("sensitive")),
                           db_label=rec.get("db_label") or f"{len(hits)} saved hit(s)",
                           hits_for_fasta="\n".join(lines),
                           demo_db=seqlab.DEMO_PROTEIN_DB)
        if rec.get("kind") == "msa":
            return _render(request, "msa.html", results=rec, error=None,
                           demo=seqlab.DEMO_MSA_SEQUENCES, chunk=60)
    return RedirectResponse(url="/history")


# ---------------------------------------------------------------------------
# Sequence laboratory: BLAST + MSA
# ---------------------------------------------------------------------------

def _parse_seq_input(text: str, upload_text: Optional[str],
                     parse_upload_fasta: bool = False) -> list[dict]:
    """Combine pasted text and uploaded FASTA into a sequence list."""
    entries = seqlab.parse_fasta(text or "")
    if upload_text:
        uploaded = seqlab.parse_fasta(upload_text)
        # Merge (avoid ids that could collide with pasted ones).
        existing = {e["id"] for e in entries}
        for e in uploaded:
            if e["id"] in existing:
                e["id"] = f"{e['id']}_{uploaded.index(e)}"
            entries.append(e)
    return entries


@dapp.get("/blast", response_class=HTMLResponse)
def blast_page(request: Request):
    return _render(request, "blast.html", results=None, error=None,
                   demo_db=seqlab.DEMO_PROTEIN_DB)


@dapp.post("/blast", response_class=HTMLResponse)
async def blast_submit(
    request: Request,
    query: str = Form(""),
    db_text: str = Form(""),
    db_type: str = Form("protein"),
    query_file: Optional[UploadFile] = File(None),
    db_file: Optional[UploadFile] = File(None),
):
    q_upload = None
    db_upload = None
    try:
        if query_file is not None:
            q_upload = (await query_file.read()).decode("utf-8", errors="replace")
        if db_file is not None:
            db_upload = (await db_file.read()).decode("utf-8", errors="replace")
    except Exception as exc:  # noqa: BLE001
        return _render(request, "blast.html", results=None, error=str(exc),
                       demo_db=seqlab.DEMO_PROTEIN_DB)

    query_seq = q_upload if q_upload else query
    query_entries = seqlab.parse_fasta(query_seq or "")
    if not query_entries:
        return _render(request, "blast.html", results=None,
                       error="Please provide a query sequence (paste it or upload a FASTA file).",
                       demo_db=seqlab.DEMO_PROTEIN_DB)
    query_str = query_entries[0]["seq"]

    if db_upload:
        database = seqlab.parse_fasta(db_upload)
    elif db_text.strip():
        database = seqlab.parse_fasta(db_text)
    else:
        # No database chosen: automatically use the built-in reference set
        # that matches the query alphabet (protein or nucleotide).
        database = (seqlab.DEMO_NUCLEOTIDE_DB if seqlab.is_nucleotide(query_str)
                    else seqlab.DEMO_PROTEIN_DB)
    db_label = seqlab.database_label(database) if hasattr(seqlab, "database_label") else f"{len(database)} sequence(s)"

    try:
        hits = seqlab.blast_search(query_str, database)
        sensitive = False
        if not hits:
            # No shared-word seeds: fall back to scoring every subject so the
            # page always returns a result table (closest matches, labelled).
            hits = seqlab.blast_search(query_str, database, require_seed=False)
            sensitive = True
    except Exception as exc:  # noqa: BLE001
        return _render(request, "blast.html", results=None,
                       error=f"Search failed: {exc}",
                       demo_db=seqlab.DEMO_PROTEIN_DB)

    rendered = [h.to_dict() for h in hits]
    db_by_id = {d["id"]: d["seq"] for d in database}

    top = hits[:10]
    lines = [f">query\n{query_str}"]
    for h in top:
        seq = db_by_id.get(h.subject_id, "")
        if seq:
            lines.append(f">{h.subject_id}\n{seq}")

    _append_seqlab({
        "kind": "blast",
        "ts": _seq_key(),
        "query": query_str,
        "db_type": "auto",
        "db_size": len(database),
        "db_label": db_label,
        "sensitive": sensitive,
        "hits": [
            {**h.to_dict(), "subject_seq": db_by_id.get(h.subject_id, "")}
            for h in hits[:50]
        ],
    })

    return _render(request, "blast.html", results=rendered, error=None,
                   query=query_str, db_size=len(database), sensitive=sensitive,
                   db_label=db_label,
                   hits_for_fasta="\n".join(lines),
                   demo_db=seqlab.DEMO_PROTEIN_DB)


@dapp.get("/msa", response_class=HTMLResponse)
def msa_page(request: Request):
    return _render(request, "msa.html", results=None, error=None,
                   demo=seqlab.DEMO_MSA_SEQUENCES)


@dapp.post("/msa", response_class=HTMLResponse)
async def msa_submit(
    request: Request,
    sequences: str = Form(""),
    seq_file: Optional[UploadFile] = File(None),
    chunk: int = Form(60),
):
    upload = None
    try:
        if seq_file is not None:
            upload = (await seq_file.read()).decode("utf-8", errors="replace")
    except Exception as exc:  # noqa: BLE001
        return _render(request, "msa.html", results=None, error=str(exc),
                       demo=seqlab.DEMO_MSA_SEQUENCES)

    text = upload if upload else sequences
    parsed = seqlab.parse_fasta(text or "")
    if len(parsed) < 2:
        return _render(request, "msa.html", results=None,
                       error="Please provide at least two sequences (FASTA).",
                       demo=seqlab.DEMO_MSA_SEQUENCES)

    try:
        result = seqlab.align_multiple(parsed)
    except Exception as exc:  # noqa: BLE001
        return _render(request, "msa.html", results=None,
                       error=f"Alignment failed: {exc}",
                       demo=seqlab.DEMO_MSA_SEQUENCES)

    data = result.to_dict()
    _append_seqlab({
        "kind": "msa",
        "ts": _seq_key(),
        "columns": data["columns"],
        "method": data["method"],
        "guide_tree": data["guide_tree"],
        "conservation": data["conservation"],
        "sequences": data["sequences"],
        "originals": [{"id": e["id"], "seq": e["seq"]} for e in parsed],
    })
    return _render(request, "msa.html", results=data, error=None,
                   demo=seqlab.DEMO_MSA_SEQUENCES, chunk=chunk)


# JSON API (useful for programmatic access).
@dapp.post("/api/blast")
async def api_blast(query: str = Form(...), db_text: str = Form(""),
                    db_type: str = Form("protein")):
    query_entries = seqlab.parse_fasta(query or "")
    if not query_entries:
        return {"error": "query is required"}
    database = (seqlab.parse_fasta(db_text) if db_text.strip()
                else (seqlab.DEMO_NUCLEOTIDE_DB if db_type == "nucleotide"
                      else seqlab.DEMO_PROTEIN_DB))
    hits = seqlab.blast_search(query_entries[0]["seq"], database)
    return {"hits": [h.to_dict() for h in hits]}


@dapp.post("/api/msa")
async def api_msa(sequences: str = Form(...)):
    parsed = seqlab.parse_fasta(sequences or "")
    if len(parsed) < 2:
        return {"error": "at least two sequences required"}
    return seqlab.align_multiple(parsed).to_dict()