"""Translate EF1% into the number of actives actually recovered in the top 1%.

EF is a ratio, so on small active sets a "doubling" can be a single compound.
This makes that explicit before any claim is written into the manuscript.
Usage: python ef_to_hits.py [autodock|diffdock]
"""
import math
import sys
from pathlib import Path

import pandas as pd

R = Path("/home/joe/projects/drug_repurposing/multipose_experiment/results")
TARGETS = ["TP53", "ESR1_ant", "PPARG"]


def main(pathway="autodock"):
    sfx = "_autodock" if pathway == "autodock" else ""
    print(f"{'target':9s} {'N':>5s} {'A':>4s} {'k':>4s}  {'arm':22s} {'EF1%':>6s} {'hits':>5s}")
    for t in TARGETS:
        b = R / f"{t}_baseline{sfx}.csv"
        if not b.exists():
            continue
        bd = pd.read_csv(b)
        N = int(bd.n.iloc[0]); A = int(bd.n_actives.iloc[0])
        k = max(1, math.ceil(0.01 * N)); prev = A / N
        hits = lambda ef: round(ef * prev * k)  # noqa: E731
        row = bd[bd.scorer == "GNINA_CNNaff"]
        if len(row):
            ef = float(row["EF1%"].iloc[0])
            print(f"{t:9s} {N:5d} {A:4d} {k:4d}  {'GNINA single-pose':22s} {ef:6.2f} {hits(ef):5d}")
        for f in [R / f"{t}_gnina_multipose{sfx}.csv", R / f"{t}_nmdn_multipose{sfx}.csv"]:
            if not f.exists():
                continue
            for _, r in pd.read_csv(f).iterrows():
                print(f"{'':9s} {'':5s} {'':4s} {'':4s}  {r['arm']:22s} "
                      f"{r['EF1%']:6.2f} {hits(r['EF1%']):5d}")
        nm = bd[bd.scorer == "NMDN"]
        if len(nm):
            ef = float(nm["EF1%"].iloc[0])
            print(f"{'':9s} {'':5s} {'':4s} {'':4s}  {'NMDN single-pose':22s} {ef:6.2f} {hits(ef):5d}")
        print()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "autodock")
