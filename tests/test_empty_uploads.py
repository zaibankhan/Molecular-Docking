"""Guard rails for empty uploaded files.

An uploaded file that arrives as 0 bytes (e.g. a local placeholder file that
was never downloaded) must produce a clear PrepError/ReceptorError instead of
RDKit's raw "OSError: File error: Invalid input file" crash.
"""
from __future__ import annotations

from pathlib import Path

from pipeline import ligand_prep as lp, receptor_prep as rp


def test_empty_ligand_sdf_is_clear_error(tmp_path: Path):
    empty = tmp_path / "empty.sdf"
    empty.write_bytes(b"")
    try:
        lp.prepare_ligand(str(empty), tmp_path / "lig.sdf", tmp_path / "lig.pdbqt", ph=7.4)
    except lp.PrepError as exc:
        assert "empty" in str(exc).lower()
    else:
        raise AssertionError("expected PrepError for empty ligand file")


def test_empty_receptor_pdb_is_clear_error(tmp_path: Path):
    empty = tmp_path / "empty.pdb"
    empty.write_bytes(b"")
    try:
        rp.clean_pdb_to_file(str(empty), tmp_path / "rec.pdb", tmp_path / "rec.pdbqt", ph=7.4)
    except rp.ReceptorError as exc:
        assert "empty" in str(exc).lower()
    else:
        raise AssertionError("expected ReceptorError for empty receptor file")


def test_smi_file_resolves_to_smiles_text(tmp_path: Path):
    smi = tmp_path / "lig.smi"
    smi.write_text("CCO\n", encoding="utf-8")
    path, smiles = lp.resolve_input(str(smi))
    assert path is None
    assert smiles == "CCO"