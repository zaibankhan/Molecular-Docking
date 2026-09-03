"""Shared input validation used by the interactive CLI and the web form.

Every user-facing input (receptor, ligand, box geometry) is validated here so the
CLI prompts and the browser UI behave consistently and give helpful feedback.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from .utils import looks_like_pdb_id

# Small helper for human-friendly parsing of comma-separated triples.
_TRIPLE_RE = re.compile(r"^\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*$")


class InputError(ValueError):
    """Raised when a user-supplied input is invalid; message is user-facing."""


def parse_triple(text: str, what: str = "box coordinate") -> Optional[tuple[float, float, float]]:
    """Parse 'X,Y,Z' into a float triple. Returns None when blank/empty."""
    text = (text or "").strip()
    if not text:
        return None
    m = _TRIPLE_RE.match(text)
    if not m:
        raise InputError(
            f"Invalid {what}: expected three comma-separated numbers like '20,20,20', got '{text}'."
        )
    return tuple(float(g) for g in m.groups())


def validate_receptor(value: str) -> None:
    """Confirm a receptor input is usable: a readable file or a plausible PDB ID."""
    value = (value or "").strip()
    if not value:
        raise InputError("Receptor is required. Give a PDB file path or a 4-character PDB ID.")
    p = Path(value)
    if p.exists():
        # Pre-empt clear mistakes.
        if p.is_dir():
            raise InputError(f"'{value}' is a folder, not a structure file. Give a .pdb file path.")
        if p.suffix.lower() not in (".pdb", ".ent", ".pdbqt"):
            # Allowed to proceed; receptor prep will decide. Keep strict about obvious issues.
            pass
    elif not looks_like_pdb_id(value):
        raise InputError(
            f"'{value}' is neither an existing file nor a 4-character RCSB PDB ID. "
            "Example: 'demo/receptor.pdb' or '1STP'."
        )


def validate_ligand(value: str) -> None:
    """Confirm a ligand input is usable: a file, a valid SMILES, or a resolvable name."""
    value = (value or "").strip()
    if not value:
        raise InputError("Ligand is required. Give a SMILES string, a file path, or a compound name.")
    p = Path(value)
    if p.exists():
        if p.is_dir():
            raise InputError(f"'{value}' is a folder. Give an .sdf/.mol/.pdb ligand file.")
        return
    # Valid SMILES (RDKit often accepts it) short-circuits before name lookup.
    from rdkit import Chem

    if Chem.MolFromSmiles(value):
        return
    # Not a file and not valid SMILES; assume a compound name (PubChem lookup).
    # We do not hard-fail here: the orchestrator will attempt the name lookup.


def validate_box(center_text: Optional[str], source_text: Optional[str]) -> None:
    """Validate optional box inputs: either an explicit center or a reference file."""
    center = parse_triple(center_text, "box center")
    if center is not None and source_text:
        raise InputError("Provide EITHER a box center OR a box-source reference file, not both.")
    if source_text:
        source_text = source_text.strip()
        if not Path(source_text).exists():
            raise InputError(
                f"Box-source file not found: '{source_text}'. Give an existing .sdf/.pdb/.pdbqt file."
            )