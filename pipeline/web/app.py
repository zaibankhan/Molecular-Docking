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
    receptor: str = Form(...),
    ligand: str = Form(...),
    exhaustiveness: int = Form(16),
    num_modes: int = Form(5),
    seed: int = Form(42),
    ph: float = Form(7.4),
    box_center: str = Form(""),
    box_source: str = Form(""),
    box_size: str = Form("20,20,20"),
):
    out_dir = None  # default -> ./runs/<timestamp>
    box_center_parts = [p.strip() for p in box_center.split(",")] if box_center.strip() else []
    if box_center_parts and len(box_center_parts) != 3:
        raise HTTPException(400, "box_center must be X,Y,Z")
    box_center_tuple = (
        tuple(float(p) for p in box_center_parts) if box_center_parts else None
    )
    box_size_tuple = tuple(float(p.strip()) for p in box_size.split(","))

    cfg = DockingConfig(
        receptor=receptor,
        ligand=ligand,
        out_dir=None,  # auto timestamped run dir
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        ph=ph,
        box_center=box_center_tuple,
        box_size=box_size_tuple,
        box_source=box_source.strip() or None,
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