"""FastAPI web interface for the molecular docking pipeline.

Run with:  python -m pipeline serve  --reload
Then open  http://127.0.0.1:8000  in your browser.
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, Request, UploadFile, HTTPException
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles

from ..config import DockingConfig


_HERE = Path(__file__).resolve().parent
_TEMPLATES_DIR = _HERE / "templates"

dapp = FastAPI(title="Molecular Docking Pipeline", version="1.0.0")
templates = Jinja2Templates(directory=str(_TEMPLATES_DIR))


@dapp.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    """Render validation/errors through the results template for a friendly UI."""
    return templates.TemplateResponse(
        request,
        "results.html",
        {"error": f"({exc.status_code}) {exc.detail}", "results": None, "out_dir": None},
        status_code=exc.status_code,
    )

# --- shared state: last finished run (so results survive page refresh) ---
_STATE: dict = {"results": None, "out_dir": None, "error": None}
_LOCK = threading.Lock()


def _latest_run_dir(base: Path) -> Optional[Path]:
    if not base.is_dir():
        return None
    candidates = sorted(base.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    for c in candidates:
        if c.is_dir() and (c / "results.json").exists():
            return c
    return None


@dapp.get("/", response_class=HTMLResponse)
def index(request: Request):
    with _LOCK:
        has_result = _STATE["results"] is not None
    return templates.TemplateResponse(
        request,
        "index.html",
        {"has_result": has_result},
    )


def _run_in_thread(cfg: DockingConfig):
    """Execute the pipeline and stash the outcome. Errors are caught so the
    server stays up even when a docking run fails."""
    from ..orchestrator import run

    try:
        results = run(cfg, verbose=False)
        out_dir = cfg.resolved_out()
        with _LOCK:
            _STATE["results"] = results
            _STATE["out_dir"] = str(out_dir)
            _STATE["error"] = None
    except Exception as exc:  # noqa: BLE001
        with _LOCK:
            _STATE["results"] = None
            _STATE["out_dir"] = None
            _STATE["error"] = f"{type(exc).__name__}: {exc}"


@dapp.post("/run")
async def do_run(
    receptor: str = Form(None),
    ligand: str = Form(None),
    box_center: str = Form(""),
    box_source: str = Form(""),
    box_size: str = Form("20,20,20"),
    exhaustiveness: int = Form(32),
    num_modes: int = Form(9),
    seed: int = Form(42),
    ph: float = Form(7.4),
    receptor_file: Optional[UploadFile] = File(None),
    ligand_file: Optional[UploadFile] = File(None),
    box_file: Optional[UploadFile] = File(None),
):
    from .. import validate as vld
    from ..utils import make_run_dir

    # Persist uploaded files into a fresh staging dir so the pipeline can read them.
    staging = make_run_dir()
    staging.mkdir(parents=True, exist_ok=True)

    async def save_upload(upload: Optional[UploadFile], name: str, allowed: tuple) -> Optional[str]:
        if upload is None:
            return None
        ext = Path(upload.filename).suffix.lower()
        if ext not in allowed:
            raise HTTPException(400, f"Unsupported file type for {name}: {ext}. Allowed: {', '.join(allowed)}")
        dest = staging / f"{name}{ext or '.txt'}"
        data = await upload.read()
        dest.write_bytes(data)
        return str(dest)

    # Resolve receptor: uploaded file wins over text field.
    receptor_path = await save_upload(
        receptor_file, "receptor", (".pdb", ".ent", ".pdbqt")
    )
    if receptor_path is None:
        try:
            vld.validate_receptor(receptor or "")
        except vld.InputError as exc:
            raise HTTPException(400, str(exc))
        receptor_value = receptor or ""
    else:
        receptor_value = receptor_path

    # Resolve ligand: uploaded file wins, else text (SMILES/path/name).
    ligand_path = await save_upload(
        ligand_file, "ligand", (".sdf", ".mol", ".pdb", ".pdbqt")
    )
    if ligand_path is None:
        try:
            vld.validate_ligand(ligand or "")
        except vld.InputError as exc:
            raise HTTPException(400, str(exc))
        ligand_value = ligand or ""
    else:
        ligand_value = ligand_path

    # Box: uploaded reference file, else text box_source / box_center.
    box_path = await save_upload(box_file, "box", (".sdf", ".pdb", ".pdbqt"))
    box_source_value = box_path if box_path else (box_source.strip() or None)
    if box_center.strip() and box_source_value:
        raise HTTPException(400, "Provide EITHER a box center OR a box-source file, not both.")
    try:
        box_center_tuple = vld.parse_triple(box_center, "box center")
    except vld.InputError as exc:
        raise HTTPException(400, str(exc))
    try:
        box_size_tuple = vld.parse_triple(box_size, "box size")
    except vld.InputError as exc:
        raise HTTPException(400, str(exc))
    if box_size_tuple is None:
        box_size_tuple = (20.0, 20.0, 20.0)

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

    thread = threading.Thread(target=_run_in_thread, args=(cfg,), daemon=True)
    thread.start()
    thread.join()

    with _LOCK:
        error = _STATE["error"]
        out_dir = _STATE["out_dir"]
    if error:
        raise HTTPException(500, error)
    return RedirectResponse(url=f"/results?out={out_dir}", status_code=303)


@dapp.get("/results", response_class=HTMLResponse)
def results(request: Request, out: Optional[str] = None):
    import json as _json

    with _LOCK:
        error = _STATE["error"]
        in_mem_results = _STATE["results"]
        in_mem_out = _STATE["out_dir"]

    # Resolve the results dir: query param > in-memory > latest on disk.
    out_dir = out or in_mem_out
    results_dict = in_mem_results
    if results_dict is None and out_dir:
        results_json = Path(out_dir) / "results.json"
        if results_json.exists():
            try:
                results_dict = _json.loads(results_json.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                results_dict = None

    if error and results_dict is None:
        return templates.TemplateResponse(
            request, "results.html", {"error": error, "results": None, "out_dir": out}
        )
    return templates.TemplateResponse(
        request,
        "results.html",
        {"results": results_dict, "out_dir": out_dir, "error": None},
    )


@dapp.get("/api/health")
def health():
    return {"status": "ok"}