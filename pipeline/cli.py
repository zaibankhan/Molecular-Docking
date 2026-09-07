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

    guide_p = sub.add_parser(
        "guided",
        help="Interactive walk-through: answer prompts to run a docking.",
    )
    guide_p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")

    doc_p = sub.add_parser("doctor", help="Check the environment and report detected tooling.")
    doc_p.add_argument("-v", "--verbose", action="store_true")

    serve_p = sub.add_parser(
        "serve",
        help="Launch a local browser-based (FastAPI) web interface.",
    )
    serve_p.add_argument("--host", default="127.0.0.1", help="Bind address (default 127.0.0.1).")
    serve_p.add_argument("--port", type=int, default=8000, help="Port (default 8000).")
    serve_p.add_argument("--reload", action="store_true", help="Auto-reload on code changes.")

    blast_p = sub.add_parser("blast", help="Local BLAST-style sequence search (fully offline).")
    blast_p.add_argument("--query", required=True, help="Query sequence (or FASTA).")
    blast_p.add_argument("--db", default="", help="Database sequences in FASTA format. "
                                                   "Empty uses the built-in demo database.")
    blast_p.add_argument("--db-type", choices=["protein", "nucleotide"], default="protein")
    blast_p.add_argument("--top", type=int, default=10, help="Max hits to report (default 10).")

    msa_p = sub.add_parser("msa", help="Multiple sequence alignment (progressive, offline).")
    msa_p.add_argument("--sequences", required=True, help="Two or more FASTA sequences.")
    msa_p.add_argument("--chunk", type=int, default=60, help="Columns per printed block (default 60).")

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


def _prompt(text: str, default: str = "") -> str:
    """Prompt for input with an optional inline default."""
    if default:
        return input(f"{text} [{default}]: ").strip() or default
    return input(f"{text}: ").strip()


def _prompt_int(text: str, default: int, minimum: int = 1) -> int:
    while True:
        raw = _prompt(text, str(default))
        try:
            val = int(raw)
        except ValueError:
            print(f"  Please enter a whole number (min {minimum}).")
            continue
        if val < minimum:
            print(f"  Value must be at least {minimum}.")
            continue
        return val


def _cmd_guided(args) -> int:
    from . import validate as vld
    from .utils import configure_logging

    print("=" * 60)
    print("Molecular Docking Pipeline - guided setup")
    print("Answer each prompt; press Enter to accept shown defaults.")
    print("=" * 60)

    # Receptor
    receptor = _prompt("Receptor (file path or 4-char PDB ID, e.g. '1STP')", "demo/receptor.pdb")
    while True:
        try:
            vld.validate_receptor(receptor)
            break
        except vld.InputError as exc:
            print(f"  {exc}")
            receptor = _prompt("Receptor (file path or 4-char PDB ID)")

    # Ligand
    ligand = _prompt(
        "Ligand (SMILES, file path, or compound name)",
        "C1C2C(C(=O)CCSC(=O)NCCC1)SCCC2C(=O)O",
    )
    while True:
        try:
            vld.validate_ligand(ligand)
            break
        except vld.InputError as exc:
            print(f"  {exc}")
            ligand = _prompt("Ligand (SMILES, file path, or compound name)")

    # Box
    box_center_text = _prompt("Box center X,Y,Z (blank to use box source or ligand centroid)", "")
    box_source = _prompt("Box-source file path (blank to skip)", "")
    while True:
        try:
            vld.validate_box(box_center_text, box_source)
            break
        except vld.InputError as exc:
            print(f"  {exc}")
            box_center_text = _prompt("Box center X,Y,Z", "")
            box_source = _prompt("Box-source file path", "")
    box_center = vld.parse_triple(box_center_text, "box center")
    box_source = box_source.strip() or None

    box_size_text = _prompt("Box size X,Y,Z (angstrom)", "20,20,20")
    while True:
        try:
            size = vld.parse_triple(box_size_text, "box size")
            if size is None:
                raise vld.InputError("Box size is required, e.g. '20,20,20'.")
            break
        except vld.InputError as exc:
            print(f"  {exc}")
            box_size_text = _prompt("Box size X,Y,Z", "20,20,20")

    # Scoring parameters
    exhaustiveness = _prompt_int("Exhaustiveness", 32)
    num_modes = _prompt_int("Number of poses to generate", 9)
    seed = _prompt_int("Random seed", 42)
    ph = 7.4
    try:
        ph_raw = _prompt("Preparation pH", "7.4")
        ph = float(ph_raw)
    except ValueError:
        ph = 7.4

    print()
    print("Starting docking with the chosen settings...")
    configure_logging(args.verbose)

    cfg = DockingConfig(
        receptor=receptor,
        ligand=ligand,
        out_dir=None,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        ph=ph,
        box_center=box_center,
        box_size=size,
        box_source=box_source,
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

    # Localhost-only by design: this web UI is intended to run on the user's own
    # device. Refuse to bind to a public/network interface to avoid accidentally
    # exposing local file operations.
    _loopback_hosts = {"127.0.0.1", "localhost", "::1"}
    host = args.host
    if host not in _loopback_hosts:
        raise SystemExit(
            f"Refusing to bind to {host!r}. The web interface is localhost-only; "
            "use 127.0.0.1 (this is intentional for security)."
        )
    port = args.port
    reload = args.reload

    url = f"http://{host}:{port}"
    print("=" * 60)
    print("Molecular Docking Pipeline - web interface")
    print(f"  Local address:  {url}")
    print(f"  Runs only on this device (localhost).")
    print(f"  Docs (OpenAPI): {url}/docs")
    print("  Press Ctrl+C to stop.")
    print("=" * 60)
    uvicorn.run("pipeline.web.app:dapp", host=host, port=port, reload=reload)
    return 0


def _cmd_blast(args) -> int:
    from . import seqlab
    from .utils import configure_logging

    configure_logging(False)
    query_entries = seqlab.parse_fasta(args.query)
    if not query_entries:
        print("Error: no query sequence found.")
        return 1
    query_str = query_entries[0]["seq"]
    if args.db.strip():
        database = seqlab.parse_fasta(args.db)
    elif args.db_type == "nucleotide":
        database = seqlab.DEMO_NUCLEOTIDE_DB
    else:
        database = seqlab.DEMO_PROTEIN_DB

    hits = seqlab.blast_search(query_str, database)[: args.top]
    print(f"BLAST-like search: {len(hits)} hit(s) against {len(database)} subject(s)")
    print(f"  query: {query_str[:40]}")
    print()
    for h in hits:
        flag = " *** significant" if h.e_value < 1e-4 else ""
        print(f"[{h.rank:>2}] {h.subject_id:<14} E={h.e_value:<10.2g} "
              f"ident={h.identity_pct:>5.1f}%  cov={h.coverage:>5.1f}%  bit={h.bit_score:>6.1f}{flag}")
        print(f"      {h.description}")
    return 0


def _cmd_msa(args) -> int:
    from . import seqlab
    from .utils import configure_logging

    configure_logging(False)
    parsed = seqlab.parse_fasta(args.sequences)
    if len(parsed) < 2:
        print("Error: provide at least two FASTA sequences.")
        return 1
    result = seqlab.align_multiple(parsed)
    print(f"Multiple sequence alignment: {len(result.sequences)} sequences, {result.columns} columns")
    print(f"Method: {result.method}")
    if result.guide_tree:
        print(f"Guide tree: {result.guide_tree}")
    print()
    for s in result.sequences:
        seq = s["seq"]
        for start in range(0, result.columns, args.chunk):
            block = seq[start : start + args.chunk]
            header = f"{s['id']:<16} {block}"
            print(header)
        print()
    if result.conservation:
        cons = "".join(
            "*" if c >= 80 else ("+" if c >= 50 else " ") for c in result.conservation
        )
        for start in range(0, result.columns, args.chunk):
            print(f"{'':<16} {cons[start:start + args.chunk]}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        return _cmd_run(args)
    if args.command == "guided":
        return _cmd_guided(args)
    if args.command == "doctor":
        return _cmd_doctor(args)
    if args.command == "serve":
        return _cmd_serve(args)
    if args.command == "blast":
        return _cmd_blast(args)
    if args.command == "msa":
        return _cmd_msa(args)
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
