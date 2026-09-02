"""Receptor preparation: PDB -> cleaned PDB + Vina-compatible PDBQT.

Strategy (dependency-light, reproducible):
  1. Parse the PDB with RDKit.
  2. Strip crystal waters and non-protein hetero residues (keeps standard
     amino acids; any small-molecule cofactor is removed by default since the
     docking box covers the protein).
  3. Add polar hydrogens at the target pH.
  4. Assign Gasteiger partial charges.
  5. Write a PDBQT with AutoDock Vina (AD4) atom types.

No MGLTools/ADFR dependency is required.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

from .utils import fetch_pdb, looks_like_pdb_id, LOG

# Standard amino acid 3-letter codes we keep in the receptor.
_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
}


class ReceptorError(RuntimeError):
    pass


def resolve_receptor_source(receptor: str, cwd: Optional[Path] = None) -> tuple[Optional[Path], Optional[str]]:
    """Resolve receptor input to (file_path_or_None, pdb_id_or_None)."""
    cwd = cwd or Path.cwd()
    as_path = Path(receptor)
    if as_path.exists():
        return as_path, None
    if looks_like_pdb_id(receptor):
        return None, receptor.strip()
    return as_path, None  # will error on read


def _is_water(resname: str) -> bool:
    return resname in ("HOH", "WAT", "H2O")


def _is_standard_aa(resname: str) -> bool:
    return resname in _AA


def clean_receptor_mol(mol: Chem.Mol) -> Chem.Mol:
    """Remove water and non-standard residues; return the protein heavy-atom mol.

    Uses per-atom PDBResidueInfo (residue groups are not exposed via
    GetResidues in current RDKit versions).
    """
    residues: dict[tuple[str, int], str] = {}
    for atom in mol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        if ri is None:
            continue
        key = (ri.GetChainId(), ri.GetResidueNumber())
        residues[key] = ri.GetResidueName().strip()

    remove_idx = []
    for atom in mol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        if ri is None:
            continue
        resname = ri.GetResidueName().strip()
        key = (ri.GetChainId(), ri.GetResidueNumber())
        # Standard amino acid residue names are kept; everything else removed.
        if resname not in _AA:
            remove_idx.append(atom.GetIdx())
    if not remove_idx:
        return mol
    rmmol = Chem.RWMol(mol)
    for idx in sorted(remove_idx, reverse=True):
        rmmol.RemoveAtom(idx)
    cleaned = rmmol.GetMol()
    try:
        Chem.SanitizeMol(cleaned)
    except Exception:
        pass
    LOG.info("Removed %d water/hetero atoms.", len(remove_idx))
    return cleaned


def add_polar_hydrogens(mol: Chem.Mol, ph: float = 7.4) -> Chem.Mol:
    """Add hydrogens; at neutral pH this is a reasonable default for proteins."""
    return Chem.AddHs(mol, addCoords=True)


# AutoDock 4 / Vina atom-type inference for common protein elements.
def _ad4_atom_type(atom: Chem.Atom) -> str:
    """Return an AutoDock (AD4) atom type for a heavy protein atom."""
    symbol = atom.GetSymbol()
    if symbol == "C":
        return "C"
    if symbol == "N":
        return "N"
    if symbol == "O":
        return "OA"  # oxygen in Vina is typed via hybridization; OA acceptable for protein O
    if symbol == "S":
        # Disulfide S vs polar S; default to SA.
        return "SA"
    if symbol == "P":
        return "P"
    if symbol == "F":
        return "F"
    if symbol == "Cl":
        return "Cl"
    if symbol == "Br":
        return "Br"
    if symbol == "I":
        return "I"
    return symbol if len(symbol) <= 2 else "C"


def clean_pdb_to_file(
    receptor: str,
    out_clean_pdb: Path,
    out_pdbqt: Path,
    ph: float = 7.4,
    cwd: Optional[Path] = None,
) -> tuple[Chem.Mol, dict]:
    """Clean a receptor and write cleaned PDB + PDBQT. Returns (mol, meta)."""
    file_path, pdb_id = resolve_receptor_source(receptor, cwd)

    if file_path is not None:
        mol = Chem.MolFromPDBFile(str(file_path), removeHs=False, sanitize=True)
        if mol is None:
            raise ReceptorError(f"Could not parse receptor PDB: {file_path}")
        source_desc = str(file_path)
    elif pdb_id is not None:
        tmp = out_clean_pdb.parent / f"_{pdb_id}_raw.pdb"
        fetch_pdb(pdb_id, tmp)
        mol = Chem.MolFromPDBFile(str(tmp), removeHs=False, sanitize=True)
        if mol is None:
            raise ReceptorError(f"Could not parse fetched PDB for {pdb_id}")
        source_desc = f"PDB:{pdb_id}"
    else:
        raise ReceptorError(f"Receptor not found and not a valid PDB ID: {receptor}")

    n_in = mol.GetNumAtoms()
    try:
        cleaned = clean_receptor_mol(mol)
    except Exception as exc:  # pragma: no cover
        raise ReceptorError(f"Receptor cleaning failed: {exc}") from exc

    cleaned = add_polar_hydrogens(cleaned, ph=ph)
    LOG.info("Receptor atoms: %d -> %d (with hydrogens).", n_in, cleaned.GetNumAtoms())

    # Assign Gasteiger charges.
    try:
        AllChem.ComputeGasteigerCharges(cleaned)
    except Exception as exc:  # pragma: no cover
        LOG.warning("Gasteiger charge assignment failed (%s); using 0.0.", exc)

    # Write cleaned PDB.
    writer = Chem.PDBWriter(str(out_clean_pdb))
    writer.write(cleaned)
    writer.close()
    LOG.info("Wrote cleaned receptor PDB: %s", out_clean_pdb)

    # Write PDBQT with AD4 atom types + Gasteiger charges + polar H only.
    _write_receptor_pdbqt(cleaned, out_pdbqt)

    aa_residues = set()
    for atom in cleaned.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        if ri is not None and _is_standard_aa(ri.GetResidueName().strip()):
            aa_residues.add((ri.GetChainId(), ri.GetResidueNumber()))
    meta = {
        "source": source_desc,
        "atoms_original": n_in,
        "atoms_prepared": cleaned.GetNumAtoms(),
        "residues": len(aa_residues),
        "ph": ph,
    }
    return cleaned, meta


def _padded(segments: list[tuple[str, int, str]]) -> str:
    """Join fixed-width PDB columns from (text, width, 'l'|'r') segments."""
    out = []
    for text, width, align in segments:
        if align == "r":
            out.append(text[:width].rjust(width))
        else:
            out.append(text[:width].ljust(width))
    return "".join(out)


def _write_receptor_pdbqt(mol: Chem.Mol, out_pdbqt: Path) -> None:
    conf = mol.GetConformer()
    lines = []
    lines.append("REMARK  Generated by molecular-docking-pipeline receptor prep")
    serial = 0
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() == "H" and not _is_polar_h(atom, mol):
            continue  # nonpolar hydrogens not needed for rigid receptor in Vina
        serial += 1
        pos = conf.GetAtomPosition(i)
        charge = atom.GetDoubleProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else 0.0
        atype = _ad4_atom_type(atom)
        pdbinfo = atom.GetPDBResidueInfo()
        resname = (pdbinfo.GetResidueName().strip() if pdbinfo else "UNK")
        resnum = (pdbinfo.GetResidueNumber() if pdbinfo else i)
        chain = (pdbinfo.GetChainId() if pdbinfo else " ")
        atomname = (pdbinfo.GetName().strip() if pdbinfo else atom.GetSymbol())
        # PDBQT column map (1-based columns):
        #   rec 1-6 ; serial 7-11 ; name 13-16 ; resname 18-21 ; chain 22 ;
        #   resseq 23-26 ; x 31-38 ; y 39-46 ; z 47-54 ; occ 55-60 ;
        #   temp 61-66 ; charge 71-76 ; type 79-80
        rec = _padded(
            [
                ("ATOM", 6, "l"),        # 1-6
                (str(serial), 5, "r"),   # 7-11 (right-aligned number)
                (" ", 1, "l"),           # 12
                (atomname, 4, "l"),      # 13-16
                (" ", 1, "l"),           # 17 altLoc
                (resname, 4, "l"),       # 18-21
                (chain or " ", 1, "l"),  # 22
                (str(resnum), 4, "r"),   # 23-26 (right-aligned number)
                (" ", 1, "l"),           # 27 iCode
                (" ", 3, "l"),           # 28-30
                (f"{float(pos.x):8.3f}", 8, "l"),  # 31-38
                (f"{float(pos.y):8.3f}", 8, "l"),  # 39-46
                (f"{float(pos.z):8.3f}", 8, "l"),  # 47-54
                ("1.00", 6, "l"),        # 55-60 occupancy
                ("0.00", 6, "l"),        # 61-66 temp factor
                (" ", 4, "l"),           # 67-70 footnote
                (f"{charge:6.2f}", 6, "l"),  # 71-76 partial charge
                (" ", 2, "l"),           # 77-78
            ]
        ) + (" " + atype)[-2:]
        lines.append(rec)
    lines.append("END")
    text = "\n".join(lines) + "\n"
    out_pdbqt.write_text(text, encoding="utf-8")
    LOG.info("Wrote receptor PDBQT: %s", out_pdbqt)


def _is_polar_h(atom: Chem.Mol, mol: Chem.Mol) -> bool:
    """Polar H = hydrogen bonded to N, O, S, or P."""
    if atom.GetSymbol() != "H":
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() in ("N", "O", "S", "P", "F", "Cl"):
            return True
    return False
