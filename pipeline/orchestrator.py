"""Pipeline orchestration: bind config -> run each stage -> emit report."""
from __future__ import annotations

from pathlib import Path

from .config import DockingConfig, ResolvedFiles
from .deps import detect
from .utils import LOG, configure_logging
from . import box as box_mod
from . import ligand_prep as lp
from . import receptor_prep as rp
from . import docking as dock
from . import analysis as an
from . import report as rep


def run(cfg: DockingConfig, verbose: bool = False) -> dict:
    """Execute the full docking pipeline for the given config."""
    configure_logging(verbose)

    env = detect()
    LOG.info("Environment:\n  %s", "\n  ".join(env.summary_lines()))
    if not env.vina_ok:
        raise RuntimeError(
            "AutoDock Vina executable was not found. Place a vina binary in the "
            "'bin/' folder or add it to PATH."
        )

    out_dir = cfg.resolved_out()
    out_dir.mkdir(parents=True, exist_ok=True)
    files = ResolvedFiles.build(out_dir)

    # --- Stage 1: ligand prep (needed early for centroid when no explicit box) ---
    LOG.info("=== Stage 1/5: Ligand preparation ===")
    ligand_meta = lp.prepare_ligand(
        cfg.ligand, files.ligand_sdf, files.ligand_pdbqt, ph=cfg.ph
    )
    ligand_mol = lp._read_mol_for_box(files.ligand_sdf)

    # --- Stage 2: receptor prep ---
    LOG.info("=== Stage 2/5: Receptor preparation ===")
    receptor_mol, receptor_meta = rp.clean_pdb_to_file(
        cfg.receptor, files.receptor_clean_pdb, files.receptor_pdbqt, ph=cfg.ph
    )

    # --- Stage 3: search box ---
    LOG.info("=== Stage 3/5: Search box ===")
    box = box_mod.default_box(
        ligand_mol,
        explicit_center=cfg.box_center,
        explicit_size=cfg.box_size,
        center_source=cfg.box_source,
    )

    # --- Stage 4: docking ---
    LOG.info("=== Stage 4/5: Docking with AutoDock Vina ===")
    docking_result = dock.run_vina(
        vina_path=env.vina_path,
        receptor_pdbqt=files.receptor_pdbqt,
        ligand_pdbqt=files.ligand_pdbqt,
        out_poses_pdbqt=files.poses_pdbqt,
        box=box,
        exhaustiveness=cfg.exhaustiveness,
        num_modes=cfg.num_modes,
        seed=cfg.seed,
        config_txt=files.config_txt,
    )

    # --- Stage 5: analysis + report ---
    LOG.info("=== Stage 5/5: Analysis & report ===")
    pose_paths = rep.split_poses_to_pdb(files.poses_pdbqt, files.pose_dir, cfg.num_modes)

    receptor_atoms = an.parse_receptor_pdbqt(files.receptor_pdbqt)
    ligand_atoms = an.parse_pose_pdbqt(files.poses_pdbqt, pose_index=1)
    interactions = an.map_interactions(receptor_atoms, ligand_atoms)

    admet = an.admet_lite(ligand_meta.get("canonical_smiles", ""))
    ligand_meta["admet"] = admet

    results = rep.build_results_dict(
        cfg, docking_result, interactions, receptor_meta, ligand_meta, box.to_dict(), pose_paths
    )
    rep.write_json(files.results_json, results)
    rep.make_report(cfg, files, docking_result, interactions, receptor_meta, ligand_meta, box.to_dict(), None)

    LOG.info("=== Done. Output in: %s ===", out_dir)
    return results
