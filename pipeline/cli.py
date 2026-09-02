"""Command-line interface for the molecular docking pipeline."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import DockingConfig
from .deps import detect, print_summary
from .utils import configure_logging


def _add_box_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--box-center",
        default=None,
        metavar="X,Y,Z",
        help="Explicit box center as three comma-separated floats (default: ligand centroid).",
    )
    parser.add_argument(
        "--box-from",
        default=None,
        metavar="FILE",
        help="Center the box on the centroid of an SDF/PDB/PDBQT file "
             "(e.g. a cocrystallized ligand or pocket-defining reference).",
    )
    parser.add_argument(
        "--box-size",
        default="20,20,20",
        metavar="X,Y,Z",
        help="Box size in angstrom as three comma-separated floats (default 20x20x20).",
    )


def _parse_triple(s: str, default=None):
    if s is None:
        return default
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"Expected 3 comma-separated values, got: {s}")
    try:
        return tuple(float(p) for p in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid numbers: {s}") from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="molecular-docking-pipeline",
        description="Receptor + ligand -> AutoDock Vina docking -> interaction report.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", help="Execute a full docking run.")
    run_p.add_argument(
        "--receptor", required=True,
        help="Receptor PDB file path OR a 4-char PDB ID (fetched from RCSB).",
    )
    run_p.add_argument(
        "--ligand", required=True,
        help="Ligand: SMILES, SDF/MOL/PDB file path, or a compound name (PubChem lookup).",
    )
    run_p.add_argument("--out", default=None, help="Output directory (default: ./runs/<timestamp>).")
    run_p.add_argument("--exhaustiveness", type=int, default=32, help="Vina exhaustiveness (default 32).")
    run_p.add_argument("--num-modes", type=int, default=9, help="Number of poses to generate (default 9).")
    run_p.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (default 42).")
    run_p.add_argument("--ph", type=float, default=7.4, help="Preparation pH (default 7.4).")
    _add_box_arg(run_p)
    run_p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")

    doc_p = sub.add_parser("doctor", help="Check the environment and report detected tooling.")
    doc_p.add_argument("-v", "--verbose", action="store_true")

    serve_p = sub.add_parser(
        "serve",
        help="Launch a local browser-based (FastAPI) web interface.",
    )
    serve_p.add_argument("--host", default="127.0.0.1", help="Bind address (default 127.0.0.1).")
    serve_p.add_argument("--port", type=int, default=8000, help="Port (default 8000).")
    serve_p.add_argument("--reload", action="store_true", help="Auto-reload on code changes.")

    return parser


def _cmd_run(args) -> int:
    box_center = _parse_triple(args.box_center)
    box_size = _parse_triple(args.box_size)
    out = Path(args.out) if args.out else None

    cfg = DockingConfig(
        receptor=args.receptor,
        ligand=args.ligand,
        out_dir=out,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        seed=args.seed,
        ph=args.ph,
        box_center=box_center,
        box_size=box_size,
        box_source=args.box_from,
    )
    from .orchestrator import run
    results = run(cfg, verbose=args.verbose)
    best = results.get("docking", {}).get("best_pose")
    print("\nBest pose affinity:", best["affinity_kcal_mol"], "kcal/mol" if best else "n/a")
    print("Output:", cfg.resolved_out())
    return 0


def _cmd_doctor(args) -> int:
    configure_logging(args.verbose)
    env = detect()
    print("Detected environment:")
    print_summary(env)
    if not env.vina_ok:
        print("\nWARNING: AutoDock Vina not found. Docking runs will fail.")
        return 1
    print("\nOK: Vina is available; pipeline is ready to run.")
    return 0


def _cmd_serve(args) -> int:
    import uvicorn  # local import so CLI stays usable without the web extras

    # Run from the project root so ./runs and relative demo/default paths resolve.
    import os

    os.chdir(str(Path(__file__).resolve().parent.parent))

    host = args.host
    port = args.port
    reload = args.reload

    url = f"http://{host}:{port}"
    print("=" * 60)
    print("Molecular Docking Pipeline — web interface")
    print(f"  Local address:  {url}")
    print(f"  Docs (OpenAPI): {url}/docs")
    print("  Press Ctrl+C to stop.")
    print("=" * 60)
    uvicorn.run("pipeline.web.app:dapp", host=host, port=port, reload=reload)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        return _cmd_run(args)
    if args.command == "doctor":
        return _cmd_doctor(args)
    if args.command == "serve":
        return _cmd_serve(args)
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
