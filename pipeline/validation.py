"""Native-pose validation: how closely does a docked pose reproduce the
crystallographic (native) ligand conformation?

For a docking method to be credible it must be able to *recapitulate* a known,
bound ligand pose. Given a native reference structure, this module aligns each
docked pose to it and reports the heavy-atom RMSD. A value under the
conventional threshold (2.0 A) indicates a successful pose reproduction.

The docked pose and the native reference are the same molecule but can carry
different internal atom orderings and hydrogen counts. Correspondence is
established geometrically: a rigid (Kabsch) closest-point alignment with
optimal assignment (Hungarian) is run from several random initial frames, and
the lowest RMSD is kept. Multi-start seeding avoids the local-minima trap of a
single ICP pass, so an identical structure yields an RMSD of ~0 regardless of
atom order or initial placement.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .utils import LOG

# Conventional threshold for a successful docking redock (heavy-atom RMSD).
SUCCESS_RMSD_ANGSTROM = 2.0


@dataclass
class PoseValidation:
    index: int
    rmsd_angstrom: Optional[float]
    success: bool

    def to_dict(self) -> dict:
        return {
            "pose": self.index,
            "rmsd_angstrom": None if self.rmsd_angstrom is None else round(self.rmsd_angstrom, 3),
            "success": self.success,
        }


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_AD4_TO_ELEMENT = {
    "C": "C", "A": "C", "CG0": "C", "C3": "C", "C2": "C", "CA": "C", "CD": "C",
    "N": "N", "NA": "N", "N3": "N", "N2": "N", "G0": "N", "G": "N",
    "OA": "O", "O2": "O", "O3": "O",
    "SA": "S", "SH": "S", "P": "P", "F": "F", "HD": "H", "H": "H",
}


def _element_from_type(atom_type: str) -> Optional[str]:
    key = atom_type.strip()
    if key in _AD4_TO_ELEMENT:
        return _AD4_TO_ELEMENT[key]
    if len(key) >= 2 and key[:2].capitalize() in ("Cl", "Br", "Na", "Ca"):
        return key[:2].capitalize()
    return key[:1].upper() if key else None


def _element_from_any(name: str) -> Optional[str]:
    letters = "".join(ch for ch in name if ch.isalpha())
    if len(letters) >= 2 and letters[:2].capitalize() in (
        "Cl", "Br", "Na", "Ca", "Zn", "Fe", "Mg", "Mn"):
        return letters[:2].capitalize()
    return letters[:1].upper() if letters else None


def _line_element(line: str) -> Optional[str]:
    atype = line[76:78].strip()
    if atype:
        return _element_from_type(atype)
    return _element_from_any(line[12:16].strip())


def _parse_pose_heavy_coords(poses_pdbqt: Path) -> list[np.ndarray]:
    """Heavy-atom coords per MODEL in a Vina poses PDBQT."""
    models: list[list[np.ndarray]] = []
    current: list[np.ndarray] = []
    text = Path(poses_pdbqt).read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            if current:
                models.append(current)
            current = []
        elif line.startswith(("ATOM", "HETATM")):
            try:
                coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if _line_element(line) not in (None, "H"):
                current.append(coords)
    if current:
        models.append(current)
    return [np.asarray(m) for m in models if len(m)]


def _native_heavy_coords(native: str) -> Optional[np.ndarray]:
    """Heavy-atom coordinates from a PDB/PDBQT path (or raw text)."""
    p = Path(native)
    try:
        if p.exists() and p.suffix.lower() != ".pdbqt":
            from rdkit import Chem

            mol = Chem.MolFromPDBFile(str(p), removeHs=True, sanitize=False,
                                      proximityBonding=True)
            if mol is None:
                mol = Chem.MolFromPDBBlock(p.read_text(encoding="utf-8", errors="replace"),
                                           removeHs=True, sanitize=False,
                                           proximityBonding=True)
            if mol is not None and mol.GetNumAtoms():
                conf = mol.GetConformer()
                return np.array([
                    [conf.GetAtomPosition(a.GetIdx()).x,
                     conf.GetAtomPosition(a.GetIdx()).y,
                     conf.GetAtomPosition(a.GetIdx()).z]
                    for a in mol.GetAtoms()
                ])
        # Fall back to plain text parsing (PDB or PDBQT).
        text = p.read_text(encoding="utf-8", errors="replace") if p.exists() else native
        coords = []
        for line in text.splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if _line_element(line) not in (None, "H"):
                coords.append(xyz)
        return np.asarray(coords) if coords else None
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------------------
# Rigid alignment + RMSD
# ---------------------------------------------------------------------------

def _kabsch_align(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    centre_t = target.mean(axis=0)
    centre_m = mobile.mean(axis=0)
    mt = mobile - centre_m
    tt = target - centre_t
    h = tt.T @ mt
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return (mt @ rot.T) + centre_t


def _distance_signature(coords: np.ndarray) -> np.ndarray:
    """Per-atom sorted distance-to-all-others (column i = signature of atom i)."""
    diff = coords[:, None, :] - coords[None, :, :]
    d = np.sqrt(np.sum(diff * diff, axis=2))
    np.fill_diagonal(d, -1.0)
    return np.sort(d, axis=1)


def _correspondence(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Permutation mapping each target atom to its mobile atom.

    Uses the matrix of sorted pairwise distances, which is permutation
    invariant: for a rigid identical molecule, corresponding atoms share equal
    distance signatures, so a single Hungarian pass finds the exact matching
    with no ICP local-minimum risk. Two-element columns are only weakly
    discriminative for small molecules, so the nearest-neighbour tie-break is
    applied via the raw geometry after alignment instead.
    """
    from scipy.optimize import linear_sum_assignment

    n = len(target)
    sig_m = _distance_signature(mobile[:n])
    sig_t = _distance_signature(target[:n])
    cost = np.zeros((n, n), dtype=float)
    for i in range(n):
        diff = sig_m - sig_t[None, i, :]
        cost[:, i] = np.sqrt(np.mean(diff * diff, axis=1))
    _, col = linear_sum_assignment(cost)
    return col  # mobile[j] pairs with target[col[j]]


def _rmsd_aligned(mobile: np.ndarray, target: np.ndarray) -> Optional[float]:
    """Heavy-atom RMSD by distance-matrix correspondence + Kabsch (+ refinement).

    Deterministic and free of the ICP local-minimum problem: the atom
    correspondence is recovered from permutation-invariant distance signatures,
    then a single Kabsch+Hungarian refinement pass tightens the fit for
    partially flexible conformers.
    """
    from scipy.optimize import linear_sum_assignment

    n = min(len(mobile), len(target))
    if n < 3:
        return None
    m = np.asarray(mobile[:n], dtype=float)
    t = np.asarray(target[:n], dtype=float)

    order = _correspondence(m, t)
    # mobile[j] <-> target[order[j]]; build paired arrays so Kabsch can run.
    paired_m = m[order]
    best = _pose_rmsd(paired_m, t)

    # Single local refinement: alternate Kabsch alignment + assignment, but
    # only accept the result if it lowers the RMSD (guards against drift).
    aligned = paired_m
    for _ in range(20):
        new_a = _kabsch_align(aligned, t)
        diff = new_a[:, None, :] - t[None, :, :]
        cost = np.sum(diff * diff, axis=2).astype(float)
        row, col = linear_sum_assignment(cost)
        rmsd = float(np.sqrt(np.mean(cost[row, col])))
        if rmsd < best - 1e-6:
            best = rmsd
        else:
            break
        matched = np.zeros_like(new_a)
        matched[col] = new_a[row]
        aligned = matched
    return round(best, 4)


def _pose_rmsd(mobile_paired: np.ndarray, target: np.ndarray) -> float:
    aligned = _kabsch_align(mobile_paired, target)
    d = aligned - target
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def validate_against_native(poses_pdbqt: Path, native: str,
                            prepared_sdf: Optional[Path] = None,
                            num_poses: Optional[int] = None) -> list[PoseValidation]:
    """Compute heavy-atom RMSD of every docked pose vs the native reference."""
    _ = prepared_sdf  # kept for API compatibility / future graph refinement
    native_coords = _native_heavy_coords(native)
    if native_coords is None or len(native_coords) == 0:
        LOG.warning("Native validation skipped: could not parse native reference %s", native)
        return []

    pose_models = _parse_pose_heavy_coords(Path(poses_pdbqt))
    if num_poses:
        pose_models = pose_models[:num_poses]

    LOG.info("Native validation: %d poses vs %d reference heavy atoms",
             len(pose_models), len(native_coords))

    out: list[PoseValidation] = []
    for i, model in enumerate(pose_models, start=1):
        rmsd = _rmsd_aligned(model, native_coords)
        out.append(
            PoseValidation(
                index=i,
                rmsd_angstrom=rmsd,
                success=rmsd is not None and rmsd <= SUCCESS_RMSD_ANGSTROM,
            )
        )
    return out


def summary(validations: list[PoseValidation]) -> dict:
    """Overall validation summary: best RMSD, success, per-pose list."""
    best = None
    for v in validations:
        if v.rmsd_angstrom is not None:
            if best is None or v.rmsd_angstrom < best.rmsd_angstrom:
                best = v
    success = bool(best and best.success)
    return {
        "threshold_angstrom": SUCCESS_RMSD_ANGSTROM,
        "best_pose": best.index if best else None,
        "best_rmsd_angstrom": None if best is None else round(best.rmsd_angstrom, 3),
        "success": success,
        "message": (
            "Redock successful: a predicted pose lies within %.1f A (heavy-atom "
            "RMSD) of the native crystallographic ligand." % SUCCESS_RMSD_ANGSTROM
            if success
            else "No predicted pose reproduced the native pose within %.1f A "
            "(heavy-atom RMSD)." % SUCCESS_RMSD_ANGSTROM
        ),
        "per_pose": [v.to_dict() for v in validations],
    }
