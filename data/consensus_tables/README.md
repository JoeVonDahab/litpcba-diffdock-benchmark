# Per-target docking + rescoring scores

One gzipped CSV per LIT-PCBA target (`<TARGET>_scores.csv.gz`), one row per ligand,
with scores from both screening pathways and the experimental activity label. These are
the processed results behind Table II and the ML re-ranking. Local file-path columns
from the working pipeline were removed; everything else is verbatim.

## Columns

| column | meaning |
|--------|---------|
| `drug` | ligand ID (PubChem CID) |
| `Active_diffdock`, `Active_autodock` | experimental activity label (True=active). Same ligand, one per pathway. |
| `target` | LIT-PCBA target name |
| **DiffDock pathway** | best-confidence DiffDock pose, then rescored |
| `confidence_score` | DiffDock confidence (native ranking score) |
| `CNNscore_diffdock`, `CNNaffinity_diffdock`, `CNN_VS_diffdock` | GNINA CNN pose-score / predicted affinity / virtual-screen score |
| `minimizedAffinity_diffdock`, `Affinity_diffdock` | GNINA (Vina) affinities (kcal/mol) |
| `NMDN-Score_diffdock` | NMDN mixture-density distance likelihood (`MDN_LOGSUM_DIST2_REFDIST2`); pose-quality score |
| `pKd-Score_diffdock` | NMDN predicted binding affinity (`PROP_PRED`, pKd) |
| `rank_nmdn_diffdock`, `rank_gnina_diffdock`, `rank_diffdock`, `mean_rank_diffdock` | precomputed ranks (1=best) |
| **AutoDock pathway** | best-affinity AutoDock pose, then rescored |
| `Affinity_kcal_per_mol` | AutoDock-GPU docking affinity (kcal/mol, lower=better) |
| `CNNscore_autodock`, `CNNaffinity_autodock`, `CNN_VS_autodock`, `minimizedAffinity_autodock`, `Affinity_autodock` | GNINA on the AutoDock pose |
| `NMDN-Score_autodock`, `pKd-Score_autodock` | NMDN distance-likelihood / predicted pKd |
| `rank_nmdn_autodock`, `rank_gnina_autodock`, `rank_autodock`, `mean_rank_autodock` | precomputed ranks |

## Scoring-method → column (for reproducing Table II)

| method | rank by | direction |
|--------|---------|-----------|
| DiffDock | `confidence_score` | higher = better |
| DiffDock-GNINA | `CNNaffinity_diffdock` | higher = better |
| DiffDock-NMDN | `pKd-Score_diffdock` | higher = better |
| AutoDock | `Affinity_kcal_per_mol` | **lower** = better |
| AutoDock-GNINA | `CNNaffinity_autodock` | higher = better |
| AutoDock-NMDN | `pKd-Score_autodock` | higher = better |

**NMDN usage.** The NMDN model emits two numbers: `NMDN-Score` (distance likelihood,
used as a pose-quality *filter* in the consensus strategies) and `pKd-Score` (predicted
affinity, used to *rank* ligands). Single-method NMDN rows rank by `pKd-Score`.

See `../../docs/targets.csv` for per-target ligand/active counts and which targets had
inactives subsampled to 5%. Reproduce with `python scripts/reproduce_metrics.py`.
