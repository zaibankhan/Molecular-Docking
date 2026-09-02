"""Ligand preparation: input parsing -> 3D -> PDBQT.

Uses RDKit for input parsing, neutralization, 3D embedding, energy
optimization and Gasteiger charges, then Meeko to write a Vina-compatible
PDBQT with a proper torsion tree.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .utils import fetch_pubchem_smiles, LOG

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    MEEKO_OK = True
except Exception as exc:  # pragma: no cover
    MEEKO_OK = False
    LOG.warning("Meeko unavailable (%s); PDBQT writing will be skipped.", exc)


class PrepError(RuntimeError):
    pass


def resolve_input(ligand: str, cwd: Optional[Path] = None) -> tuple[Optional[Path], str]:
    """Resolve a ligand input to (file_path_or_None, smiles).

    Resolution order: existing SDF/MOL file -> SMILES -> compound name via PubChem.
    """
    cwd = cwd or Path.cwd()
    as_path = Path(ligand)
    if as_path.exists() and as_path.suffix.lower() in (".sdf", ".mol", ".pdb"):
        return as_path, ""

    mol = Chem.MolFromSmiles(ligand.strip())
    if mol is not None:
        return None, ligand.strip()

    # Fall through to name lookup.
    LOG.info("Input is not a file or valid SMILES; treating as compound name.")
    return None, fetch_pubchem_smiles(ligand.strip())


def _neutralize(mol: Chem.Mol) -> Chem.Mol:
    """Remove common fragment charges to yield a neutral species where possible."""
    patts = (
        Chem.MolFromSmarts("[N+;H2]"),
        Chem.MolFromSmarts("[N+;H1]"),
        Chem.MolFromSmarts("[N+;H0]"),
        Chem.MolFromSmarts("[O-]"),
        Chem.MolFromSmarts("[N+]#C"),
        Chem.MolFromSmarts("[S+]"),
    )
    working = Chem.Mol(mol)
    for patt in patts:
        if patt is None:
            continue
        while working.HasSubstructMatch(patt):
            rms = AllChem.ReplaceSubstructs(
                working, patt, Chem.MolFromSmiles("[*]"), replaceAll=True
            )
            if not rms:
                break
            new = Chem.Mol(rms[0])
            Chem.SanitizeMol(new)
            working = new
    return working


def prepare_ligand(
    ligand_input: str,
    out_sdf: Path,
    out_pdbqt: Path,
    ph: float = 7.4,
    cwd: Optional[Path] = None,
) -> dict:
    """Prepare a ligand and write SDF + PDBQT. Returns a metadata dict."""
    if Chem is None:
        raise PrepError("RDKit is not installed; cannot prepare ligand.")

    file_path, smiles = resolve_input(ligand_input, cwd)

    if file_path is not None:
        if file_path.suffix.lower() == ".sdf":
            supplier = Chem.SDMolSupplier(str(file_path), removeHs=False, sanitize=True)
            mol = next((m for m in supplier if m is not None), None)
            if mol is None:
                raise PrepError(f"No valid molecule in SDF file: {file_path}")
        elif file_path.suffix.lower() == ".mol":
            mol = Chem.MolFromMolFile(str(file_path), sanitize=True, removeHs=False)
        elif file_path.suffix.lower() == ".pdb":
            mol = Chem.MolFromPDBFile(str(file_path), removeHs=False, sanitize=True)
        else:
            raise PrepError(f"Unsupported ligand file type: {file_path.suffix}")
        if mol is None:
            raise PrepError(f"Could not parse ligand file: {file_path}")
    else:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if mol is None:
            raise PrepError(f"Invalid SMILES: {smiles}")

    mol = _neutralize(mol)

    # Embed 3D coordinates.
    mol3 = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol3, params) != 0:
        # Retry with a different random seed and enforce chirality.
        params2 = AllChem.ETKDGv3()
        params2.randomSeed = 0xBEEF
        params2.enforceChirality = True
        if AllChem.EmbedMolecule(mol3, params2) != 0:
            raise PrepError("3D embedding failed for ligand.")

    # Optimize geometry.
    try:
        AllChem.MMFFOptimizeMolecule(mol3, maxIters=1000)
    except Exception as exc:  # pragma: no cover
        LOG.warning("UFF/MMFF optimization failed (%s); using embedded conformer.", exc)

    # Compute canonical (non-3D) descriptors for stability.
    canon = Chem.MolToSmiles(Chem.RemoveHs(mol))
    heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # PDBQT via Meeko (Gasteiger charges handled internally).
    if MEEKO_OK:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol3)
        if not mol_setups:
            raise PrepError("Meeko produced no torsion setups for the ligand.")
        pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            raise PrepError(f"Meeko PDBQT writing failed: {err}")
        out_pdbqt.write_text(pdbqt_string, encoding="utf-8")
        LOG.info("Wrote ligand PDBQT: %s", out_pdbqt)
    else:
        raise PrepError("Meeko is required to write the ligand PDBQT; install 'meeko'.")

    writer = Chem.SDWriter(str(out_sdf))
    writer.write(mol3)
    writer.close()
    LOG.info("Wrote ligand SDF: %s", out_sdf)

    meta = {
        "input": ligand_input,
        "source": "file" if file_path is not None else ("smiles" if not file_path and smiles != "" else "name"),
        "canonical_smiles": canon,
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "mol_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "heavy_atoms": heavy,
        "rotatable_bonds": rotatable,
        "n_atoms": mol3.GetNumAtoms(),
    }
    return meta


def _read_mol_for_box(sdf_path: Path) -> "Chem.Mol":
    """Read a prepared ligand SDF back as an RDKit mol (with coords) for box centering."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise PrepError(f"Could not re-read ligand SDF for box: {sdf_path}")
    return mol


def smiles_coordinates(smiles: str) -> tuple[float, float, float]:
    """Return the 3D centroid of a molecule from its SMILES (for box centering)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise PrepError("Could not embed molecule to compute box center.")
    conf = mol.GetConformer()
    xs, ys, zs = [], [], []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(p.x)
        ys.append(p.y)
        zs.append(p.z)
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    return (round(float(cx), 2), round(float(cy), 2), round(float(cz), 2))
