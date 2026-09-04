"""FastAPI web interface for the molecular docking pipeline.

Purpose: a clean, reliable browser front-end to the CLI/pipeline engine.
The user fills in receptor + ligand + box, clicks "Run docking", and is shown
a progress page while Vina works, then lands on a results page.

Run with:  python -m pipeline serve
Then open http://127.0.0.1:8000 in your browser. Localhost-only by design.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from ..config import DockingConfig
from .. import validate as vld

_HERE = Path(__file__).resolve().parent
_TEMPLATES_DIR = _HERE / "templates"
_RUNS_BASE = Path("runs")

dapp = FastAPI(title="Molecular Docking Pipeline", version="1.3.0")
templates = Jinja2Templates(directory=str(_TEMPLATES_DIR))

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
    """Write an uploaded file to the staging dir. Returns its path or None."""
    if upload is None:
        return None
    filename = upload.filename or ""
    ext = Path(filename).suffix.lower()
    if ext not in allowed:
        raise ValueError(
            f"Unsupported file type for {name}: '{ext or 'no-extension'}'. "
            f"Allowed: {', '.join(allowed)}"
        )
    data = await upload.read()
    if len(data) > _MAX_UPLOAD:
        raise ValueError(f"Uploaded {name} file is too large (>100 MB).")
    dest = staging / f"{name}_{staging_token()}{ext}"
    dest.write_bytes(data)
    return str(dest)


def staging_token() -> str:
    import time

    return str(int(time.time() * 1000))


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
    box_file: Optional[UploadFile] = File(None),
):
    errors: list[str] = []

    receptor = (receptor or "").strip()
    ligand = (ligand or "").strip()
    box_center = (box_center or "").strip()
    box_source = (box_source or "").strip()
    box_size = (box_size or "").strip()

    if not receptor and receptor_file is None:
        errors.append("Receptor is required (upload a PDB file or type a path / PDB ID).")
    if not ligand and ligand_file is None:
        errors.append("Ligand is required (upload a file or type SMILES / path / name).")

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

    prev = {"receptor": receptor, "ligand": ligand, "box_center": box_center,
            "box_source": box_source, "box_size": box_size,
            "exhaustiveness": exhaustiveness, "num_modes": num_modes, "seed": seed, "ph": ph}

    if errors:
        return _render(request, "index.html", running=False, errors=errors, prev=prev)

    # --- Persist uploaded files. ---
    staging = _RUNS_BASE / "_uploads"
    staging.mkdir(parents=True, exist_ok=True)
    try:
        rec_path = await _persist_upload(receptor_file, "receptor", (".pdb", ".ent", ".pdbqt"), staging)
        lig_path = await _persist_upload(ligand_file, "ligand", (".sdf", ".mol", ".pdb", ".pdbqt"), staging)
        box_path = await _persist_upload(box_file, "box", (".sdf", ".pdb", ".pdbqt"), staging)
    except ValueError as exc:
        return _render(request, "index.html", running=False, errors=[str(exc)], prev=prev)

    receptor_value = rec_path or receptor
    ligand_value = lig_path or ligand
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
    return _render(request, "running.html", running=True)


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