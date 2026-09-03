"""Tests for the shared input-validation layer."""
from __future__ import annotations

import pytest

from pipeline import validate as vld


def test_receptor_valid_path():
    vld.validate_receptor("demo/receptor.pdb")  # should not raise


def test_receptor_valid_pdb_id():
    vld.validate_receptor("1STP")


def test_receptor_invalid():
    with pytest.raises(vld.InputError):
        vld.validate_receptor("a path that does not exist and is not an id!")


def test_ligand_smiles():
    vld.validate_ligand("OC(=O)c1ccccc1C(=O)O")


def test_ligand_name_passes_as_plausible():
    vld.validate_ligand("aspirin")  # treated as a name; no hard-fail


def test_ligand_fails_only_when_truly_empty():
    with pytest.raises(vld.InputError):
        vld.validate_ligand("   ")


def test_parse_triple_ok():
    assert vld.parse_triple("20,20,20") == (20.0, 20.0, 20.0)
    assert vld.parse_triple(" 1.5 , -2 , 3 ") == (1.5, -2.0, 3.0)


def test_parse_triple_bad():
    with pytest.raises(vld.InputError):
        vld.parse_triple("a,b,c")


def test_box_center_and_source_conflict():
    with pytest.raises(vld.InputError):
        vld.validate_box("1,2,3", "demo/pocket_btn.pdb")


def test_box_source_missing_file():
    with pytest.raises(vld.InputError):
        vld.validate_box("", "does/not/exist.sdf")