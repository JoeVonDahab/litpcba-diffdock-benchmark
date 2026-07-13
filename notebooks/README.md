# Notebooks — consensus and ML reproduction (author's original code)

These are the analysis notebooks used for the paper, with outputs cleared. They
reproduce the parts of Table II / Table III that depend on run-specific choices
(rank weights, NMDN-score threshold), which is why they are shipped as code
rather than re-derived in `scripts/`.

## Inputs

Both notebooks read the per-target merged tables. In the deposit these live at
`../data/consensus_tables/<TARGET>__merged.csv.gz` (global intersection of both
pathways) and `<TARGET>__autodock_consensus.csv.gz` (AutoDock-pathway consensus).
Point the `path1…path15` / merged-file variables at those files
(`pd.read_csv(..., compression="gzip")`).

## `consensus_reproduction.ipynb`

Rank-based consensus (CC-Medium / UC-Strong / CC-Weak, per-pathway and GLOBAL).
Core function `merge_and_rank`: inner-join NMDN ∩ GNINA ∩ dock on `drug`, rank
each method by its score (NMDN→`pKd-Score`, GNINA→`CNNaffinity`, dock→
`confidence_score`, higher = better), then sort by the weighted mean rank.
`rank_globally` averages the six ranks across both pathways. EF/BEDROC/ROC-AUC
are computed on the resulting ranking.

## `ml_reranking_table3.ipynb` — Table III

Supervised re-ranking (WNN, XGBoost, LightGBM/LambdaMART, RandomForest, DeepMLP)
on the merged feature tables. Needs an env with `lightgbm`, `xgboost`, `optuna`,
`torch` (the project's `what` conda env).

### Important caveat (carry into the manuscript)

Reported ML EF (Table III) is computed on the **validation subset only**, under a
`GroupShuffleSplit` on `subgroup_id`, where each large target is chopped into
≤10k-row subgroups. Subgroups of the **same target** therefore land in both train
and validation → this is a **within-target** split, **not** held-out-target. The
models see target-specific labels. Table III (ML, val subset, same targets) is
**not** on the same footing as Table II (docking/consensus, full library) and
must not be presented as a fair head-to-head.
