# Per-method score tables (Table II single-method source)

One gzipped CSV per (target, method): `<TARGET>__<METHOD>.csv.gz`, where METHOD is
`DiffDock`, `DiffDock-GNINA`, `DiffDock-NMDN`, `AutoDock-GNINA`, or `AutoDock-NMDN`.
Each file is that method's **full ligand set** for the target — the authoritative
source for the single-method rows of Table II / Supplementary S2.

## Columns

| column | meaning |
|--------|---------|
| `ligand_id` | ligand ID (PubChem CID) |
| `score` | the method's native ranking score, higher = better (DiffDock→confidence, GNINA→CNNaffinity, NMDN→pKd-Score) |
| `active` | experimental activity label (1 = active, 0 = inactive) |

## Reproduce

```python
import pandas as pd
from scripts.metrics import ef_at
d = pd.read_csv("data/per_target_scores/OPRK1__AutoDock-GNINA.csv.gz")
ef_at(d.active.values, d.score.values, 0.01)   # -> 12.5
```

`scripts/reproduce_table2.py` does this across all targets/methods. The merged
tables (with every raw score column, for consensus + ML) are in
`../consensus_tables/`.
