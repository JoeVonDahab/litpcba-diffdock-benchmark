"""Parse GNINA-scored multi-pose SDFs -> per-ligand aggregated scores -> EF.
Arms (GNINA CNNaffinity across the 20 poses):
  B_max     : max CNNaffinity over poses   (standard multi-pose "best score")
  C_top3    : mean of top-3 CNNaffinity     (robust variant)
  A ref     : the single-pose baseline (from baseline.py)
Usage: python aggregate_gnina.py TP53
"""
import sys, re
import numpy as np, pandas as pd
from pathlib import Path
from metrics import all_metrics

BENCH = Path("/home/joe/projects/drug_repurposing/jobs/benchmarking")
EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")

CNNAFF = re.compile(r"<CNNaffinity>\s*\n\s*([-\d.eE]+)")


def parse_scored(sdf_text):
    """Return list of CNNaffinity floats, one per model in the scored SDF."""
    return [float(x) for x in CNNAFF.findall(sdf_text)]


def main(target):
    scored_dir = EXP / "work" / target / "scored"
    files = sorted(scored_dir.glob("*_scored.sdf"))
    if not files:
        print(f"no scored SDFs in {scored_dir} — run run_gnina_multipose.sh first")
        sys.exit(1)

    rec = {}  # drugid -> array of CNNaffinity per pose
    for f in files:
        drugid = f.name.replace("_scored.sdf", "")
        aff = parse_scored(f.read_text())
        if aff:
            rec[drugid] = np.array(aff, float)

    # labels
    lab = pd.read_csv(BENCH / target / f"merged_{target}.csv")
    lab["drug"] = lab["drug"].astype(str)
    label = dict(zip(lab["drug"], lab["Active_diffdock"].astype(bool).astype(int)))

    drugs = [d for d in rec if d in label]
    y = np.array([label[d] for d in drugs])
    b_max = np.array([rec[d].max() for d in drugs])
    c_top3 = np.array([np.sort(rec[d])[::-1][:3].mean() for d in drugs])

    print(f"{target}: {len(drugs)} ligands scored, {y.sum()} actives "
          f"(avg {np.mean([len(rec[d]) for d in drugs]):.1f} poses/ligand)\n")
    rows = []
    for arm, s in [("B_max_over20", b_max), ("C_mean_top3", c_top3)]:
        m = all_metrics(y, s, higher_better=True)
        rows.append({"arm": arm, "scorer": "GNINA_CNNaff", **m})
        print(f"  {arm:14s} EF1%={m['EF1%']:.3f}  EF10%={m['EF10%']:.3f}  "
              f"AUC={m['ROC_AUC']:.3f}  BEDROC={m['BEDROC20']:.4f}")
    out = EXP / "results" / f"{target}_gnina_multipose.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53")
