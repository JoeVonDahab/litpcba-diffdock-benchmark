"""Export the consensus + ML input tables into the deposit, sanitized (drop local
file-path columns) and gzipped. These feed:
  - consensus reproduction (CC-Medium/UC-Strong/CC-Weak) via notebooks/final_analysis
  - ML re-ranking / Table III via notebooks/new_algorithm

Source: jobs/benchmarking/<T>/merged_<T>.csv  (global intersection, both pathways)
        jobs/benchmarking/<T>/[*_autodock_output/]for_merge_autodock.csv (AutoDock consensus)
"""
import pandas as pd, glob, os
HERE = os.path.dirname(__file__)
BENCH = "/home/joe/projects/drug_repurposing/jobs/benchmarking"
OUT = os.path.join(HERE, "..", "data", "consensus_tables")
TARGETS = ["ADRB2","ALDH1","ESR1_ago","ESR1_ant","FEN1","GBA","IDH1","KAT2A",
           "MAPK1","MTORC1","OPRK1","PKM2","PPARG","TP53","VDR"]
DROP = ("ligand", "ligand_diffdock", "ligand_autodock", "original_file",
        "DLG_file", "smiles", "SMILES")


def sanitize(df):
    return df[[c for c in df.columns if c not in DROP]]


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in glob.glob(os.path.join(OUT, "*.csv.gz")):
        os.remove(f)
    for t in TARGETS:
        for src, tag in [(f"{BENCH}/{t}/merged_{t}.csv", "merged"),
                         (glob.glob(f"{BENCH}/{t}/for_merge_autodock.csv") +
                          glob.glob(f"{BENCH}/{t}/*_autodock_output/for_merge_autodock.csv"),
                          "autodock_consensus")]:
            path = src if isinstance(src, str) else (src[0] if src else None)
            if not path or not os.path.exists(path):
                print(f"  skip {t} {tag} (missing)")
                continue
            df = sanitize(pd.read_csv(path, low_memory=False))
            dest = os.path.join(OUT, f"{t}__{tag}.csv.gz")
            df.to_csv(dest, index=False, compression="gzip")
            print(f"  {t:9s} {tag:18s} rows={len(df):>7d} cols={df.shape[1]}")


if __name__ == "__main__":
    main()
