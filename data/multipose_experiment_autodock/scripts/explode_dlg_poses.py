"""Explode AutoDock DLG runs into per-pose SDFs for multi-pose rescoring.

AutoDock-GPU writes ONE .dlg per ligand holding all N independent runs (10 here).
The published protocol keeps only the lowest-energy run; this keeps all of them so
GNINA/NMDN can be run in their intended multi-pose mode on the AutoDock pathway.

Output: work_autodock/<T>/pose_sdfs/<drugid>__pose<k>.sdf   (k=1 is best energy)
Usage:  python explode_dlg_poses.py TP53 [n_poses]

ponytail: one obabel call per ligand (-m splits models) instead of one per pose.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/home/joe/projects/drug_repurposing")
BENCH = ROOT / "jobs" / "benchmarking"
EXP = ROOT / "multipose_experiment"

# Resolve obabel without needing the conda env activated.
OBABEL = shutil.which("obabel") or str(Path.home() / "miniconda3/envs/diffdock_nmdn/bin/obabel")

MODEL_RE = re.compile(r"DOCKED: MODEL\s+(\d+)\s*\n(.*?)\nDOCKED: ENDMDL", re.DOTALL)
ENERGY_RE = re.compile(r"Estimated Free Energy of Binding\s*=\s*([-+]?\d*\.?\d+)")


def dlg_dir(target):
    """The one *_autodock_output/.../gpu*_results dir holding this target's DLGs."""
    hits = sorted(BENCH.glob(f"{target}/*_autodock_output/autodock_dlgs_output/*/"))
    hits = [h for h in hits if any(h.glob("*.dlg"))]
    if not hits:
        raise FileNotFoundError(f"no DLG directory with *.dlg for {target}")
    return hits[0]


def parse_dlg(path, n_poses):
    """Return up to n_poses (energy, clean_pdbqt) sorted best (lowest) energy first."""
    text = path.read_text(errors="ignore")
    poses = []
    for _model, body in MODEL_RE.findall(text):
        m = ENERGY_RE.search(body)
        if not m:
            continue
        clean = "\n".join(
            l.replace("DOCKED: ", "")
            for l in body.split("\n")
            if l.strip() and not l.startswith("DOCKED: USER")
        )
        poses.append((float(m.group(1)), clean))
    poses.sort(key=lambda p: p[0])
    return poses[:n_poses]


def explode_one(args):
    dlg, outdir, n_poses = args
    drug = Path(dlg).stem
    try:
        poses = parse_dlg(Path(dlg), n_poses)
        if not poses:
            return drug, 0, "no poses"
        # Skip work already done (resumable, like the DiffDock arm).
        if (Path(outdir) / f"{drug}__pose1.sdf").exists():
            return drug, len(poses), "cached"
        with tempfile.TemporaryDirectory() as td:
            pq = Path(td) / f"{drug}.pdbqt"
            with open(pq, "w") as fh:
                for i, (energy, body) in enumerate(poses, 1):
                    fh.write(f"MODEL {i}\nREMARK binding energy {energy:.2f} kcal/mol\n")
                    fh.write(body.rstrip("\n") + "\nENDMDL\n")
            stem = Path(td) / f"{drug}__pose.sdf"
            r = subprocess.run(
                [OBABEL, str(pq), "-O", str(stem), "-osdf", "-m"],
                capture_output=True, text=True,
            )
            made = sorted(Path(td).glob(f"{drug}__pose*.sdf"))
            if not made:
                return drug, 0, f"obabel failed: {r.stderr.strip()[:120]}"
            for f in made:
                if f.stat().st_size:
                    os.replace(f, Path(outdir) / f.name)
        return drug, len(made), "ok"
    except Exception as e:  # one bad ligand must not kill the sweep
        return drug, 0, f"error: {e}"


def main(target, n_poses=10):
    src = dlg_dir(target)
    outdir = EXP / "work_autodock" / target / "pose_sdfs"
    outdir.mkdir(parents=True, exist_ok=True)
    dlgs = sorted(str(p) for p in src.glob("*.dlg"))
    print(f"[{target}] {len(dlgs)} DLGs in {src} -> up to {n_poses} poses each")
    n_lig = n_pose = 0
    bad = []
    with ProcessPoolExecutor(max_workers=min(32, (os.cpu_count() or 4))) as pool:
        futs = [pool.submit(explode_one, (d, str(outdir), n_poses)) for d in dlgs]
        for i, fut in enumerate(as_completed(futs), 1):
            drug, k, status = fut.result()
            if k:
                n_lig += 1
                n_pose += k
            else:
                bad.append((drug, status))
            if i % 500 == 0:
                print(f"  {i}/{len(dlgs)} ligands, {n_pose} poses", flush=True)
    print(f"[{target}] {n_lig} ligands -> {n_pose} pose SDFs in {outdir}")
    if bad:
        print(f"[{target}] {len(bad)} failed, e.g. {bad[:3]}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53",
         int(sys.argv[2]) if len(sys.argv) > 2 else 10)
