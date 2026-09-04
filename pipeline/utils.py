"""Shared helpers: run dirs, network I/O, logging."""
from __future__ import annotations

import datetime as _dt
import logging
import sys
from pathlib import Path
from typing import Optional
from urllib.request import urlopen, Request

LOG = logging.getLogger("pipeline")


def _noop(*_a, **_k):
    pass


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    pipeline = logging.getLogger("pipeline")
    # Only attach our console handler once; never destroy handlers that other
    # frameworks (e.g. uvicorn) have already installed.
    if not any(getattr(h, "_pipeline_console", False) for h in pipeline.handlers):
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        handler._pipeline_console = True  # noqa: SLF001
        pipeline.addHandler(handler)
    pipeline.setLevel(level)
    pipeline.propagate = False
    # Quiet noisy third-party loggers.
    for noisy in ("urllib3", "PIL"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


def make_run_dir(base: Optional[Path] = None) -> Path:
    base = base or Path.cwd() / "runs"
    base.mkdir(parents=True, exist_ok=True)
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = base / stamp
    # Ensure uniqueness for concurrent runs.
    i = 1
    while out.exists():
        out = base / f"{stamp}_{i}"
        i += 1
    out.mkdir(parents=True)
    (out / "poses").mkdir(parents=True)
    return out


def fetch_url_text(url: str, timeout: float = 60.0) -> str:
    req = Request(url, headers={"User-Agent": "moleculardocking-pipeline/1.0"})
    with urlopen(req, timeout=timeout) as resp:  # noqa: S310 (allow-listed hosts)
        return resp.read().decode("utf-8", errors="replace")


# PDB identifier regex: 1 alphanumeric + 3 alphanumeric (upper).
PDB_ID_RE = r"^[A-Za-z0-9]{4}$"


def looks_like_pdb_id(s: str) -> bool:
    import re

    return (
        bool(re.fullmatch(PDB_ID_RE, s.strip()))
        and not s.endswith((".pdb", ".PDB"))
        and "/" not in s
        and "\\" not in s
    )


def fetch_pdb(pdb_id: str, out_path: Path) -> Path:
    """Download a PDB file from RCSB into out_path; return path."""
    pdb_id = pdb_id.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    LOG.info("Fetching PDB %s from RCSB...", pdb_id)
    text = fetch_url_text(url)
    if "HEADER" not in text and "ATOM" not in text and "CRYST1" not in text:
        raise RuntimeError(f"RCSB returned no structure for PDB id '{pdb_id}'")
    out_path.write_text(text, encoding="utf-8")
    LOG.info("Saved receptor to %s", out_path)
    return out_path


def fetch_pubchem_smiles(name: str, timeout: float = 60.0) -> str:
    """Look up a compound name via PubChem PUG and return an isomeric SMILES."""
    from urllib.parse import quote

    safe = quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe}/property/CanonicalSMILES/JSON"
    LOG.info("Looking up compound '%s' on PubChem...", name)
    import json

    text = fetch_url_text(url, timeout=timeout)
    payload = json.loads(text)
    try:
        smiles = payload["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except (KeyError, IndexError, TypeError) as exc:
        raise RuntimeError(f"Could not resolve compound name '{name}' on PubChem") from exc
    return smiles
