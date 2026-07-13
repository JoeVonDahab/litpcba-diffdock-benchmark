"""Aggregate per-pose NMDN scores -> per-ligand arms -> EF.
Input: work/<T>/nmdn_multipose.csv  (lig path, NMDN-Score, pKd-Score), lig = <drugid>__pose<k>.sdf
Arms:
  NMDN_max    : max NMDN-Score over poses          (multi-pose "best score")
  NMDN_top3   : mean of top-3 NMDN-Score           (robust variant)
  pKd_bestNMDN: pKd of the NMDN-selected pose        (NMDN pose-selector -> affinity readout; Arm D)
Usage: python aggregate_nmdn.py TP53
"""
import sys, re
import numpy as np, pandas as pd
from pathlib import Path
from metrics import all_metrics

BENCH = Path("/home/joe/projects/drug_repurposing/jobs/benchmarking")
EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")


def main(target):
    csv = EXP / "work" / target / "nmdn_multipose.csv"
    if not csv.exists() or csv.stat().st_size == 0:
        print(f"no NMDN scores at {csv} — run run_nmdn_multipose.sh first")
        sys.exit(1)
    df = pd.read_csv(csv)
    df["drug"] = df["lig"].apply(lambda p: Path(p).name.split("__")[0])

    lab = pd.read_csv(BENCH / target / f"merged_{target}.csv")
    lab["drug"] = lab["drug"].astype(str)
    label = dict(zip(lab["drug"], lab["Active_diffdock"].astype(bool).astype(int)))

    g = df.groupby("drug")
    drugs, y, nmdn_max, nmdn_top3, pkd_best = [], [], [], [], []
    for d, sub in g:
        if d not in label:
            continue
        drugs.append(d); y.append(label[d])
        nm = sub["NMDN-Score"].values
        nmdn_max.append(nm.max())
        nmdn_top3.append(np.sort(nm)[::-1][:3].mean())
        pkd_best.append(sub.loc[sub["NMDN-Score"].idxmax(), "pKd-Score"])
    y = np.array(y)
    print(f"{target}: {len(drugs)} ligands, {y.sum()} actives, "
          f"avg {g.size().mean():.1f} poses/ligand\n")
    rows = []
    for arm, s in [("NMDN_max_over20", nmdn_max), ("NMDN_mean_top3", nmdn_top3),
                   ("pKd_at_bestNMDN", pkd_best)]:
        m = all_metrics(y, np.array(s), higher_better=True)
        rows.append({"arm": arm, "scorer": "NMDN", **m})
        print(f"  {arm:16s} EF1%={m['EF1%']:.3f}  EF10%={m['EF10%']:.3f}  "
              f"AUC={m['ROC_AUC']:.3f}  BEDROC={m['BEDROC20']:.4f}")
    out = EXP / "results" / f"{target}_nmdn_multipose.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53")
