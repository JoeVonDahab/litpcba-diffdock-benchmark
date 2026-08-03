"""Concatenate each ligand's per-pose AutoDock SDFs into one multi-model SDF.

NMDN reads one molecule per file (pose_sdfs/), but GNINA scores every model in a
single file, so it needs the combined form — same split the DiffDock arm uses.
Input:  work_autodock/<T>/pose_sdfs/<drugid>__pose<k>.sdf
Output: work_autodock/<T>/allposes/<drugid>.sdf   (models ordered pose1..poseN)
Usage:  python assemble_autodock_poses.py TP53
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")


def main(target):
    posedir = EXP / "work_autodock" / target / "pose_sdfs"
    outdir = EXP / "work_autodock" / target / "allposes"
    outdir.mkdir(parents=True, exist_ok=True)
    by_drug = defaultdict(list)
    for p in posedir.glob("*__pose*.sdf"):
        drug, k = p.name.split("__pose")
        by_drug[drug].append((int(k[:-4]), p))  # strip ".sdf"; env is py3.8, no removesuffix
    n = 0
    for drug, poses in by_drug.items():
        poses.sort()  # pose1 (best energy) first
        text = "".join(p.read_text() for _k, p in poses)
        if text.strip():
            (outdir / f"{drug}.sdf").write_text(text)
            n += 1
    print(f"[{target}] assembled {n} ligands -> {outdir}")
    if n:
        counts = [len(v) for v in by_drug.values()]
        print(f"[{target}] poses/ligand: min {min(counts)} max {max(counts)} "
              f"mean {sum(counts)/len(counts):.1f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53")
