"""Search-box computation.

The docking box defines the region Vina samples poses in. By default the box
is centered on the prepared ligand's centroid (blind-ish local docking) with a
configurable size. Callers may override center/size explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .utils import LOG


@dataclass
class Box:
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    description: str = "centered on ligand centroid"

    def vina_args(self) -> list[str]:
        cx, cy, cz = self.center
        sx, sy, sz = self.size
        return [
            "--center_x", f"{cx:.2f}",
            "--center_y", f"{cy:.2f}",
            "--center_z", f"{cz:.2f}",
            "--size_x", f"{sx:.2f}",
            "--size_y", f"{sy:.2f}",
            "--size_z", f"{sz:.2f}",
        ]

    def to_dict(self) -> dict:
        return {
            "center": list(self.center),
            "size": list(self.size),
            "description": self.description,
        }


def ligand_centroid(mol, ) -> tuple[float, float, float]:
    """Compute the geometric centroid of a prepared ligand molecule."""
    if mol is None:
        raise ValueError("A molecule is required to compute the box center.")
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
    LOG.info("Ligand centroid (box center): (%.2f, %.2f, %.2f)", cx, cy, cz)
    return (float(cx), float(cy), float(cz))


def centroid_from_file(path: Path) -> tuple[float, float, float]:
    """Compute the centroid of a molecule from an SDF/PDB/PDBQT file.

    Useful for centering the docking box on a known binding site (e.g. a
    cocrystallized reference ligand or a pocket-defining residue cluster).
    """
    from rdkit import Chem

    mol = None
    suffix = path.suffix.lower()
    if suffix in (".sdf", ".mol"):
        suppl = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        mol = next((m for m in suppl if m is not None), None)
    elif suffix in (".pdb", ".pdbqt", ".ent"):
        mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
        if mol is None:
            # PDBQT files are not always readable by RDKit; fall back to a
            # lightweight coordinate parser for PDBQT.
            if suffix == ".pdbqt":
                return _centroid_pdbqt(path)
    if mol is None:
        raise ValueError(f"Could not read molecule file for box centering: {path}")
    return ligand_centroid(mol)


def _centroid_pdbqt(path: Path) -> tuple[float, float, float]:
    xs, ys, zs = [], [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        raise ValueError(f"No atoms in PDBQT file: {path}")
    return (
        float(sum(xs) / len(xs)),
        float(sum(ys) / len(ys)),
        float(sum(zs) / len(zs)),
    )


def default_box(
    mol,
    explicit_center: Optional[tuple[float, float, float]] = None,
    explicit_size: Optional[tuple[float, float, float]] = None,
    center_source: Optional[str] = None,
) -> Box:
    """Compute a box. Priority: explicit_center > center_source file > ligand centroid."""
    if explicit_center is not None:
        center = explicit_center
        desc = "explicit center"
    elif center_source:
        center = centroid_from_file(Path(center_source))
        desc = f"centered on {Path(center_source).name}"
    else:
        center = ligand_centroid(mol)
        desc = "centered on ligand centroid"
    size = explicit_size if explicit_size is not None else (20.0, 20.0, 20.0)
    box = Box(center=tuple(float(v) for v in center), size=tuple(float(v) for v in size))
    box.description = desc
    return box
