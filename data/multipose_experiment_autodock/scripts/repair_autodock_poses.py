"""Valence-repair the AutoDock pose SDFs for NMDN, matching the production pipeline.

The production AutoDock pipeline converts PDBQT twice:
  convert_pdbqt_to_sdf_cli.py                     -> poses fed to GNINA   (plain obabel)
  convert_pdbqt_to_sdf_cli_with_ligand_fixing.py  -> poses fed to NMDN    (valence-repaired)

obabel's PDBQT->SDF routinely produces over-valent atoms, which RDKit refuses to
sanitize, so NMDN must be given the repaired form or ~23% of poses are discarded.
This applies the SAME repair function the production script uses.

Input:  work_autodock/<T>/pose_sdfs/        (plain obabel, left untouched for GNINA)
Output: work_autodock/<T>/pose_sdfs_fixed/  (valence-repaired, for NMDN)
Usage:  python repair_autodock_poses.py TP53
"""
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path("/home/joe/projects/drug_repurposing")
EXP = ROOT / "multipose_experiment"
sys.path.insert(0, str(ROOT / "AutoDOCK"))
sys.path.insert(0, str(ROOT / "DiffDock-NMDN"))

from convert_pdbqt_to_sdf_cli_with_ligand_fixing import repair_sdf_file  # noqa: E402


def one(args):
    src, dst = args
    try:
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)
        ok, _msg, _log = repair_sdf_file(dst)
        if not ok:
            return 0
        # Must pass the same gate nmdn_drive applies, else it is dropped anyway.
        from filter_valid_sdfs import is_valid_sdf
        return 1 if is_valid_sdf(dst) else 0
    except Exception:
        return 0


def main(target):
    src_dir = EXP / "work_autodock" / target / "pose_sdfs"
    out_dir = EXP / "work_autodock" / target / "pose_sdfs_fixed"
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(str(p), str(out_dir / p.name)) for p in src_dir.glob("*.sdf")]
    print(f"[{target}] repairing {len(jobs)} poses -> {out_dir}", flush=True)
    good = 0
    with ProcessPoolExecutor(max_workers=min(24, os.cpu_count() or 4)) as pool:
        for i, r in enumerate(pool.map(one, jobs, chunksize=200), 1):
            good += r
            if i % 10000 == 0:
                print(f"  {i}/{len(jobs)} ... {good} valid", flush=True)
    pct = 100.0 * good / len(jobs) if jobs else 0.0
    print(f"[{target}] {good}/{len(jobs)} poses valid after repair ({pct:.1f}%)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53")
