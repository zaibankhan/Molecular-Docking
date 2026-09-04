"""AutoDock Vina docking engine wrapper.

Invokes the bundled vina binary as a subprocess with a fixed config, then
parses the log and the poses file into structured results.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .utils import LOG


@dataclass
class Pose:
    index: int
    affinity: float  # kcal/mol (lower = stronger)
    rmsd_ub: Optional[float]
    rmsd_lb: Optional[float]

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "affinity_kcal_mol": round(self.affinity, 2),
            "rmsd_upper_bound": self.rmsd_ub,
            "rmsd_lower_bound": self.rmsd_lb,
        }


@dataclass
class DockingResult:
    poses: list[Pose] = field(default_factory=list)
    vina_version: str = ""
    log: str = ""
    raw_pdbqt: Path = None

    def best(self) -> Optional[Pose]:
        if not self.poses:
            return None
        return sorted(self.poses, key=lambda p: p.affinity)[0]

    def to_dict(self) -> dict:
        return {
            "num_modes": len(self.poses),
            "best_pose": self.best().to_dict() if self.best() else None,
            "poses": [p.to_dict() for p in self.poses],
            "vina_version": self.vina_version,
        }


# Matches Vina pose lines:  <mode#>  <affinity>  <rmsd l.b.>  <rmsd u.b.>
_POSE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+([0-9.]+)\s+([0-9.]+)\s*$")


def parse_log(log_text: str, expected_modes: int) -> tuple[list[Pose], str]:
    poses: list[Pose] = []
    version = ""
    for line in log_text.splitlines():
        ls = line.strip()
        if ls.lower().startswith("autodock vina"):
            version = ls
        m = _POSE_RE.match(ls)
        if m:
            try:
                idx = int(m.group(1))
                aff = float(m.group(2))
            except ValueError:
                continue
            poses.append(
                Pose(
                    index=idx,
                    affinity=aff,
                    rmsd_lb=float(m.group(3)) if m.group(3) else None,
                    rmsd_ub=float(m.group(4)) if m.group(4) else None,
                )
            )
        if len(poses) >= expected_modes:
            break
    return poses[:expected_modes], version


def run_vina(
    vina_path: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_poses_pdbqt: Path,
    box,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
    config_txt: Optional[Path] = None,
) -> DockingResult:
    """Run Vina on prepared receptor/ligand within `box`. Returns parsed result."""
    if not Path(vina_path).exists():
        raise FileNotFoundError(f"Vina executable not found: {vina_path}")

    if config_txt is not None:
        lines = [
            f"receptor = {receptor_pdbqt}",
            f"ligand = {ligand_pdbqt}",
            f"out = {out_poses_pdbqt}",
            f"center_x = {box.center[0]:.2f}",
            f"center_y = {box.center[1]:.2f}",
            f"center_z = {box.center[2]:.2f}",
            f"size_x = {box.size[0]:.2f}",
            f"size_y = {box.size[1]:.2f}",
            f"size_z = {box.size[2]:.2f}",
            f"exhaustiveness = {exhaustiveness}",
            f"num_modes = {num_modes}",
            f"seed = {seed}",
        ]
        config_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
        LOG.info("Wrote Vina config: %s", config_txt)

    cmd = [
        vina_path,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(out_poses_pdbqt),
        "--center_x", f"{box.center[0]:.2f}",
        "--center_y", f"{box.center[1]:.2f}",
        "--center_z", f"{box.center[2]:.2f}",
        "--size_x", f"{box.size[0]:.2f}",
        "--size_y", f"{box.size[1]:.2f}",
        "--size_z", f"{box.size[2]:.2f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
    ]
    LOG.info("Running Vina: %s", " ".join(cmd))
    # NOTE: stdin must be DEVNULL (not inherited). Inheriting an unusable or
    # open-but-unread stdin is a common cause of Vina hanging when launched from
    # a server / job context on Windows.
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, stdin=subprocess.DEVNULL)
    stdout = proc.stdout or ""
    stderr = proc.stderr or ""
    LOG.info("Vina exit code: %s", proc.returncode)
    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed (exit {proc.returncode}):\n{stderr or stdout}")

    full_log = (stdout + "\n" + stderr).strip()
    poses, version = parse_log(full_log, num_modes)
    if not poses:
        raise RuntimeError(
            "Vina produced no docking poses. Check receptor/ligand PDBQT validity.\n"
            + full_log[-2000:]
        )

    LOG.info("Docking complete: %d poses, best affinity %.2f kcal/mol",
             len(poses), poses[0].affinity if poses else float("nan"))
    return DockingResult(poses=poses, vina_version=version, log=full_log, raw_pdbqt=out_poses_pdbqt)
