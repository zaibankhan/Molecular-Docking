"""Pipeline configuration and run parameters."""
from __future__ import annotations

import socket
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


from .utils import make_run_dir


@dataclass
class DockingConfig:
    """Immutable-style configuration describing a single docking run."""

    # Inputs
    receptor: str                # PDB file path or 4-char PDB ID
    ligand: str                  # SMILES, SDF file path, or compound name
    out_dir: Optional[Path] = None

    # Docking engine
    exhaustiveness: int = 32
    num_modes: int = 9
    seed: int = 42
    scoring: str = "vina"

    # Search box (optional; default derived from ligand centroid)
    box_center: Optional[tuple[float, float, float]] = None
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0)
    box_source: Optional[str] = None  # SDF/PDB/PDBQT file to center the box on

    # Prep knobs
    ph: float = 7.4

    # Memoized resolved output directory (avoid creating multiple dirs).
    _out: Optional[Path] = field(default=None, repr=False)

    def resolved_out(self) -> Path:
        if self._out is None:
            self._out = self.out_dir if self.out_dir is not None else make_run_dir()
        return self._out

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("receptor", None)
        d.pop("ligand", None)
        d.pop("_out", None)
        d["out_dir"] = str(self.resolved_out())
        d["hostname"] = socket.gethostname()
        return d


@dataclass
class ResolvedFiles:
    """Concrete file paths bound to a prepared run directory."""

    out_dir: Path
    receptor_clean_pdb: Path
    receptor_pdbqt: Path
    ligand_sdf: Path
    ligand_pdbqt: Path
    config_txt: Path
    poses_pdbqt: Path
    results_json: Path
    report_md: Path
    pose_dir: Path

    @classmethod
    def build(cls, out_dir: Path) -> "ResolvedFiles":
        return cls(
            out_dir=out_dir,
            receptor_clean_pdb=out_dir / "receptor.pdb",
            receptor_pdbqt=out_dir / "receptor.pdbqt",
            ligand_sdf=out_dir / "ligand.sdf",
            ligand_pdbqt=out_dir / "ligand.pdbqt",
            config_txt=out_dir / "config.txt",
            poses_pdbqt=out_dir / "poses.pdbqt",
            results_json=out_dir / "results.json",
            report_md=out_dir / "report.md",
            pose_dir=out_dir / "poses",
        )
