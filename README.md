# Benchmarking single-pose docking, consensus rescoring, and supervised ML on LIT-PCBA

Processed data and analysis code for the study benchmarking **DiffDock, AutoDock-GPU,
GNINA, and DiffDock-NMDN** on 15 LIT-PCBA targets (manuscript AAPSJ-D-26-00302).
This repository is the reproducibility deposit: per-method docking/rescoring scores,
the consensus and ML input tables, the multi-pose control experiment, and scripts +
notebooks that regenerate every reported number.

## What's here

```
data/
  per_target_scores/     <TARGET>__<METHOD>.csv.gz — [ligand_id, score, active]
                         one file per method per target (the raw per-method full
                         ligand set). Authoritative source for Table II single-method rows.
  consensus_tables/      <TARGET>__merged.csv.gz (global intersection, both pathways;
                         also the ML feature table) and <TARGET>__autodock_consensus.csv.gz
  multipose_experiment/  single- vs multi-pose rescoring results (GNINA + NMDN)
docs/
  RECONCILIATION.md      cell-by-cell audit of Table II / Supp S2 (read this)
  supplementary_s2_corrected.csv corrected per-target EF1% values reported in S2
  table2_reproduced.csv  clean per-target EF1/EF10/ROC-AUC/BEDROC matrix (from raw tables)
  targets.csv            per-target ligand/active counts + subsampling flag
scripts/
  metrics.py                 EF1/EF10/ROC-AUC/BEDROC (rdkit-backed); has a self-check
  reproduce_table2.py        regenerate Table II single-method EF1% from the raw tables
  build_release_data.py      rebuild data/per_target_scores + docs/table2_reproduced.csv
  build_consensus_ml_data.py rebuild data/consensus_tables
notebooks/
  consensus_reproduction.ipynb  CC-Medium/UC-Strong/CC-Weak consensus (author's code)
  ml_reranking_table3.ipynb     Table III supervised re-ranking (author's code)
  README.md                     inputs + the within-target-split caveat
```

## Reproduce

```bash
pip install -r requirements.txt

# Table II — single-method rows (authoritative, independent recompute)
python scripts/reproduce_table2.py

# Full metric matrix + rebuild the gzipped score tables
python scripts/build_release_data.py     # -> docs/table2_reproduced.csv

# Consensus (Table II consensus rows) and Table III ML: run the notebooks
#   against data/consensus_tables/*.csv.gz  (see notebooks/README.md)
```

`reproduce_table2.py` and `build_release_data.py` point `BENCH` at the working
data tree; to run purely from the deposit, point them at `data/per_target_scores/`.

## Table II accuracy — audited

The single-method columns were recomputed independently and **reproduce the paper
exactly, 14/15**. The one exception is a spreadsheet transcription error in the
**ADRB2** row (its AutoDock-GNINA `9.07` and AutoDock CC-Medium `7.25` were copied
from GBA — both values are mathematically impossible for ADRB2's 17-active count).
Corrected: ADRB2 AutoDock-GNINA `9.07 → 0.00`, AutoDock CC-Medium `7.25 → 5.88`,
AutoDock-GNINA median EF1% `2.14 → 2.03`. Full proof and per-cell table in
[`docs/RECONCILIATION.md`](docs/RECONCILIATION.md).

**NMDN scores.** NMDN emits `NMDN-Score` (mixture-density distance likelihood; used
as a pose-quality filter in the consensus strategies) and `pKd-Score` (predicted
affinity; used to rank ligands in the single-method NMDN rows).

**ML caveat (Table III).** ML EF is on the validation subset under a within-target
`GroupShuffleSplit`, **not** held-out-target — not a fair head-to-head with the
full-library docking rows. See `notebooks/README.md`.

## Multi-pose control

`data/multipose_experiment/` tests whether rescoring all 20 DiffDock poses (vs one)
helps early enrichment. EF1% is unchanged in 5 of 6 target×scorer cases; the one
gain is GNINA on ESR1_ant (2.22→4.44). See its `FINDINGS.md`.

## Data source & license

Built on **LIT-PCBA** (Tran-Nguyen et al., J Chem Inf Model 2020,
doi:10.1021/acs.jcim.0c00155). Code: MIT (`LICENSE`). Derived data (`data/`): CC-BY-4.0.

## Citation

See `CITATION.cff`. Please cite the manuscript and this dataset DOI once minted.
