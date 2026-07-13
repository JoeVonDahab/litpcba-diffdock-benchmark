"""Recompute the single-method Table II / Supplementary-S2 EF metrics from the RAW
per-method output tables (each method's full ligand set), labels from actives.smi.
This is the authoritative reproduction: single-method rows use the per-method tables,
NOT the merged inner-join (which is only for consensus).

Ranking recipe (from the analysis notebooks):
  DiffDock       molecule_results_ranked.csv   confidence_score  (higher)
  *-GNINA        gnina_output.csv              CNNaffinity       (higher)
  *-NMDN         nmdn_output.csv               pKd-Score         (higher; NMDN-Score = quality filter)

Point BENCH at the working data tree. Prints per-target EF1% and medians.
"""
import pandas as pd, numpy as np, glob, os

BENCH = "/home/joe/projects/drug_repurposing/jobs/benchmarking"
TARGETS = ["ADRB2","ALDH1","ESR1_ago","ESR1_ant","FEN1","GBA","IDH1","KAT2A",
           "MAPK1","MTORC1","OPRK1","PKM2","PPARG","TP53","VDR"]


def active_ids(t):
    s = set()
    for p in glob.glob(f"{BENCH}/{t}/**/actives.smi", recursive=True):
        for line in open(p):
            w = line.split()
            if len(w) >= 2:
                s.add(w[-1].strip())
        break
    return s


def ids(df):
    for c in ("drug", "molecule_name"):
        if c in df:
            return df[c].astype(str)
    return df["lig"].astype(str).apply(lambda p: os.path.basename(str(p)).split("_best_pose")[0].split(".")[0])


def ef(df, col, acts, frac=0.01, higher=True):
    a = ids(df).isin(acts).values
    n = len(df); na = int(a.sum())
    if na == 0 or col not in df:
        return np.nan
    order = np.argsort(-df[col].values if higher else df[col].values)
    k = max(1, int(n * frac))
    return float((a[order][:k].sum() / k) / (na / n))


def autodir(t):
    d = glob.glob(f"{BENCH}/{t}/*_autodock_output")
    return d[0] if d else None


def main():
    METH = {  # method -> (file relative to pathway dir, score col)
        "GNINA": ("gnina_output.csv", "CNNaffinity"),
        "NMDN":  ("nmdn_output.csv",  "pKd-Score"),
    }
    rows = []
    for t in TARGETS:
        acts = active_ids(t)
        dd, ad = f"{BENCH}/{t}", autodir(t)
        # DiffDock baseline
        try:
            m = pd.read_csv(f"{dd}/molecule_results_ranked.csv", low_memory=False)
            m["confidence_score"] = pd.to_numeric(m["confidence_score"], errors="coerce")
            rows.append((t, "DiffDock", ef(m.dropna(subset=["confidence_score"]), "confidence_score", acts)))
        except Exception as e:
            rows.append((t, "DiffDock", np.nan))
        for path, base in [("DiffDock", dd), ("AutoDock", ad)]:
            if not base:
                continue
            for meth, (fn, col) in METH.items():
                try:
                    df = pd.read_csv(f"{base}/{fn}", low_memory=False)
                    rows.append((t, f"{path}-{meth}", ef(df, col, acts)))
                except Exception:
                    rows.append((t, f"{path}-{meth}", np.nan))
    res = pd.DataFrame(rows, columns=["target", "method", "EF1"])
    piv = res.pivot(index="target", columns="method", values="EF1")
    order = ["DiffDock", "DiffDock-GNINA", "DiffDock-NMDN", "AutoDock-GNINA", "AutoDock-NMDN"]
    piv = piv[[c for c in order if c in piv.columns]]
    print(piv.round(3).to_string())
    print("\nMedian across targets:")
    print(piv.median().round(3).to_string())
    piv.to_csv(os.path.join(os.path.dirname(__file__), "..", "docs", "table2_from_raw.csv"))


if __name__ == "__main__":
    main()
