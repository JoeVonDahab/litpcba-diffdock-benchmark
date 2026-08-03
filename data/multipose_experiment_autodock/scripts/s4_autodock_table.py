"""Build the AutoDock-pathway multi-pose control table (Supplementary Table S4).

Mirrors Table S3 (the DiffDock control) so the two are directly comparable:
per target x rescorer, the single-pose value vs the best value over the evaluated
multi-pose aggregation rules, for EF1% and EF10%.

Arm A comes from baseline.py <T> autodock; the multi-pose arms from
aggregate_nmdn.py / aggregate_gnina.py <T> autodock.

Usage: python s4_autodock_table.py
"""
import math
import sys
from pathlib import Path

import pandas as pd

EXP = Path("/home/joe/projects/drug_repurposing/multipose_experiment")
R = EXP / "results"
TARGETS = ["TP53", "ESR1_ant", "PPARG"]

# Display names for the aggregation rules, keyed by the arm labels the
# aggregators emit. Same vocabulary as Table S3.
RULES = {
    "NMDNscore_max_over10": "maximum",
    "NMDNscore_mean_top3": "mean top 3",
    "pKd_at_bestNMDN": "pKd at best NMDN",
    "pKd_max_over10": "pKd maximum",
    "pKd_mean_top3": "pKd mean top 3",
    "B_max_over10": "maximum",
    "C_mean_top3": "mean top 3",
}
# The paper ranks NMDN by NMDN-Score, so the NMDN-Score arms are the primary
# comparison; the pKd arms are the pose-selector variant (Arm D).
SCORER_ARMS = {
    "NMDN": ["NMDNscore_max_over10", "NMDNscore_mean_top3", "pKd_at_bestNMDN"],
    "GNINA": ["B_max_over10", "C_mean_top3"],
}
BASELINE_ROW = {"NMDN": "NMDN", "GNINA": "GNINA_CNNaff"}


def load(target):
    """Return (baseline_df, multi_df) for one target, or None if incomplete."""
    b = R / f"{target}_baseline_autodock.csv"
    if not b.exists():
        return None, None
    multi = []
    for f in (R / f"{target}_nmdn_multipose_autodock.csv",
              R / f"{target}_gnina_multipose_autodock.csv"):
        if f.exists():
            multi.append(pd.read_csv(f))
    return pd.read_csv(b), (pd.concat(multi, ignore_index=True) if multi else None)


def main():
    rows, missing = [], []
    for target in TARGETS:
        base, multi = load(target)
        if base is None:
            missing.append(f"{target}: no baseline")
            continue
        for scorer, arms in SCORER_ARMS.items():
            bsel = base[base.scorer == BASELINE_ROW[scorer]]
            if multi is None or bsel.empty:
                missing.append(f"{target}/{scorer}")
                continue
            msel = multi[multi.arm.isin(arms)]
            if msel.empty:
                missing.append(f"{target}/{scorer}")
                continue
            n = int(bsel["n"].iloc[0]); n_act = int(bsel["n_actives"].iloc[0])
            cells = {"target": target, "rescorer": scorer, "n": n, "n_actives": n_act}
            k = max(1, math.ceil(0.01 * n))
            for metric in ["EF1%", "EF10%"]:
                single = round(float(bsel[metric].iloc[0]), 2)
                shown = msel[metric].round(2)
                best = float(shown.max())
                cells[f"single {metric}"] = single
                cells[f"best multi {metric}"] = best
                cells[f"delta {metric}"] = round(best - single, 2)
                cells[f"rule {metric}"] = " / ".join(
                    RULES.get(a, a) for a in msel.arm[shown == best])
                if metric == "EF1%":
                    # EF is a ratio; on 24-90 actives a "doubling" can be one compound.
                    to_hits = lambda ef: round(ef * (n_act / n) * k)  # noqa: E731
                    cells["actives in top 1%"] = f"{to_hits(single)} → {to_hits(best)} of {k}"
            rows.append(cells)

    if not rows:
        print("nothing to combine yet:", "; ".join(missing))
        sys.exit(1)
    df = pd.DataFrame(rows)
    out = R / "combined_autodock.csv"
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {out}")
    if missing:
        print("INCOMPLETE — still missing: " + "; ".join(missing))


if __name__ == "__main__":
    main()
