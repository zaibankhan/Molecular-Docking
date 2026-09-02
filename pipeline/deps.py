"""Environment detection for the docking pipeline.

Locates the AutoDock Vina executable, verifies Python tooling, and reports
which components are available. Centralizes assumptions so callers can fail
early with helpful messages.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class Environment:
    python: str
    python_version: str
    rdkit: Optional[str] = None
    meeko: Optional[str] = None
    vina_path: Optional[str] = None
    vina_version: Optional[str] = None

    @property
    def vina_ok(self) -> bool:
        return bool(self.vina_path and self.vina_version)

    def summary_lines(self) -> list[str]:
        lines = [
            f"Python      : {self.python_version} ({self.python})",
            f"RDKit       : {self.rdkit or 'MISSING'}",
            f"Meeko       : {self.meeko or 'MISSING'}",
            f"AutoDock Vina: {self.vina_version or 'NOT FOUND'} ({self.vina_path or '-'})",
        ]
        return lines


def _pkg_version(name: str) -> Optional[str]:
    try:
        mod = __import__(name)
        version = getattr(mod, "__version__", None)
        return str(version) if version else "installed"
    except Exception:
        return None


def _find_vina_candidates(project_root: Path) -> list[str]:
    """Return a ranked list of candidate vina executable paths."""
    candidates: list[str] = []
    # 1. Explicit bundled binary.
    for name in ("vina.exe", "vina"):
        p = project_root / "bin" / name
        if p.exists():
            candidates.append(str(p))
    # 2. On PATH.
    found = shutil.which("vina")
    if found:
        candidates.append(found)
    # 3. Common install locations.
    common = [
        r"C:\Program Files\AutoDock Vina\vina.exe",
        r"C:\Program Files (x86)\AutoDock Vina\vina.exe",
        r"C:\AutoDockVina\vina.exe",
        r"C:\vina\vina.exe",
    ]
    for c in common:
        if Path(c).exists():
            candidates.append(c)
    # De-duplicate, preserve order.
    seen: set[str] = set()
    ordered: list[str] = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered


def _vina_version(path: str) -> Optional[str]:
    try:
        out = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        out_text = (out.stdout or "") + (out.stderr or "")
        for line in out_text.splitlines():
            if "AutoDock Vina" in line or "Vina" in line:
                return line.strip()
        return out_text.strip() or None
    except Exception:
        return None


def detect(project_root: Optional[Path] = None) -> Environment:
    project_root = project_root or (Path(__file__).resolve().parent.parent)

    env = Environment(
        python=sys.executable,
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        rdkit=_pkg_version("rdkit"),
        meeko=_pkg_version("meeko"),
    )

    for cand in _find_vina_candidates(project_root):
        ver = _vina_version(cand)
        if ver:
            env.vina_path = cand
            env.vina_version = ver
            break

    return env


def print_summary(env: Environment) -> None:
    for line in env.summary_lines():
        print(f"  {line}")
