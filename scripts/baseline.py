"""Arm A: reproduce the published SINGLE-POSE EF from merged_<TARGET>.csv.
No scoring needed — validates our metric code against the paper's numbers.
Usage: python baseline.py TP53
"""
import sys, pandas as pd
from pathlib import Path
from metrics import all_metrics

BENCH = Path("/home/joe/projects/drug_repurposing/jobs/benchmarking")
EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")

# (column, higher_is_better) for each single-pose scorer on the DiffDock pathway
DIFFDOCK_SCORERS = {
    "DiffDock_conf": ("confidence_score", True),
    "GNINA_CNNaff": ("CNNaffinity_diffdock", True),
    "NMDN": ("NMDN-Score_diffdock", True),
}


def main(target):
    df = pd.read_csv(BENCH / target / f"merged_{target}.csv")
    y = df["Active_diffdock"].astype(bool).astype(int).values
    print(f"{target}: {len(df)} ligands, {y.sum()} actives ({y.mean()*100:.2f}%)\n")
    rows = []
    for name, (col, hb) in DIFFDOCK_SCORERS.items():
        if col not in df.columns:
            print(f"  skip {name}: no column {col}")
            continue
        m = all_metrics(y, df[col].values, higher_better=hb)
        m = {"arm": "A_single_pose", "scorer": name, **m}
        rows.append(m)
        print(f"  {name:14s} EF1%={m['EF1%']:.3f}  EF10%={m['EF10%']:.3f}  "
              f"AUC={m['ROC_AUC']:.3f}  BEDROC={m['BEDROC20']:.4f}")
    out = EXP / "results" / f"{target}_baseline.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TP53")
