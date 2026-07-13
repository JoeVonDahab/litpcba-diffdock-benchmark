"""Build the deposit from the RAW per-method output tables (each method's full
ligand set per target/pathway) — the authoritative single-method source.

For every (target, method) it writes data/per_target_scores/<target>__<method>.csv.gz
with columns [ligand_id, score, active] (sanitized: id, score, 0/1 label only),
and computes full metrics (EF1%, EF10%, ROC-AUC, BEDROC20) into
docs/table2_reproduced.csv. Labels come from actives.smi.

Ranking recipe (from the analysis notebooks):
  DiffDock     molecule_results_ranked.csv  confidence_score
  *-GNINA      gnina_output.csv             CNNaffinity
  *-NMDN       nmdn_output.csv              pKd-Score   (NMDN-Score is the consensus quality filter)
"""
import pandas as pd, numpy as np, glob, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from metrics import all_metrics

HERE = os.path.dirname(__file__)
BENCH = "/home/joe/projects/drug_repurposing/jobs/benchmarking"
OUT = os.path.join(HERE, "..", "data", "per_target_scores")
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


def id_series(df):
    for c in ("drug", "molecule_name"):
        if c in df:
            return df[c].astype(str)
    return df["lig"].astype(str).apply(lambda p: os.path.basename(str(p)).split("_best_pose")[0].split(".")[0])


def autodir(t):
    d = glob.glob(f"{BENCH}/{t}/*_autodock_output")
    return d[0] if d else None


def sources(t):
    dd, ad = f"{BENCH}/{t}", autodir(t)
    src = [
        ("DiffDock",       dd, "molecule_results_ranked.csv", "confidence_score"),
        ("DiffDock-GNINA", dd, "gnina_output.csv",            "CNNaffinity"),
        ("DiffDock-NMDN",  dd, "nmdn_output.csv",             "pKd-Score"),
    ]
    if ad:
        src += [
            ("AutoDock-GNINA", ad, "gnina_output.csv", "CNNaffinity"),
            ("AutoDock-NMDN",  ad, "nmdn_output.csv",  "pKd-Score"),
        ]
    return src


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(OUT, "*.csv.gz")):
        os.remove(f)
    rows = []
    for t in TARGETS:
        acts = active_ids(t)
        for meth, base, fn, col in sources(t):
            path = os.path.join(base, fn)
            if not os.path.exists(path):
                continue
            df = pd.read_csv(path, low_memory=False)
            if col not in df.columns:
                continue
            out = pd.DataFrame({
                "ligand_id": id_series(df).values,
                "score": pd.to_numeric(df[col], errors="coerce").values,
            }).dropna(subset=["score"])
            out["active"] = out["ligand_id"].isin(acts).astype(int)
            if out["active"].sum() == 0:
                continue
            out.to_csv(os.path.join(OUT, f"{t}__{meth}.csv.gz"), index=False,
                       compression="gzip")
            m = all_metrics(out["active"].values, out["score"].values, higher_better=True)
            rows.append({"target": t, "method": meth, **m})
            print(f"{t:9s} {meth:15s} EF1%={m['EF1%']:.3f}  n={m['n']:>7d}  na={m['n_actives']}")
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(HERE, "..", "docs", "table2_reproduced.csv"), index=False)
    piv = tab.pivot(index="target", columns="method", values="EF1%")
    order = ["DiffDock","DiffDock-GNINA","DiffDock-NMDN","AutoDock-GNINA","AutoDock-NMDN"]
    piv = piv[[c for c in order if c in piv.columns]]
    print("\nEF1% median across targets:")
    print(piv.median().round(3).to_string())


if __name__ == "__main__":
    main()
