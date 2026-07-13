"""Combine baseline (Arm A) + NMDN multipose (+ GNINA multipose if present) into
one comparison table across targets. Writes results/combined.csv and prints it.
"""
import pandas as pd
from pathlib import Path

EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")
R = EXP / "results"
TARGETS = ["TP53", "ESR1_ant", "PPARG"]

frames = []
for t in TARGETS:
    for suffix, tag in [("baseline", "A_single"),
                        ("nmdn_multipose", "NMDN_multi"),
                        ("gnina_multipose", "GNINA_multi")]:
        f = R / f"{t}_{suffix}.csv"
        if f.exists():
            d = pd.read_csv(f)
            d.insert(0, "target", t)
            d.insert(1, "group", tag)
            frames.append(d)

comb = pd.concat(frames, ignore_index=True)
cols = ["target", "group", "arm", "scorer", "EF1%", "EF10%", "ROC_AUC", "BEDROC20", "n", "n_actives"]
comb = comb[[c for c in cols if c in comb.columns]]
comb.to_csv(R / "combined.csv", index=False)

# focused NMDN single-vs-multi EF1% pivot
print("\n=== EF1%: single-pose (A) vs multi-pose NMDN, per target ===")
for t in TARGETS:
    base = comb[(comb.target == t) & (comb.arm == "A_single_pose") & (comb.scorer == "NMDN")]
    mx = comb[(comb.target == t) & (comb.arm == "NMDN_max_over20")]
    if len(base) and len(mx):
        print(f"  {t:10s} single={base['EF1%'].iloc[0]:6.2f}  "
              f"multi_max={mx['EF1%'].iloc[0]:6.2f}  "
              f"Δ={mx['EF1%'].iloc[0]-base['EF1%'].iloc[0]:+.2f}")
print(f"\nwrote {R/'combined.csv'}  ({len(comb)} rows)")
