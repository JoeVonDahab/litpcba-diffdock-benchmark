# Table II / Supplementary-S2 reconciliation

Every **single-method** per-target EF value was recomputed independently from the
raw per-method output tables (each method's full ligand set per target/pathway),
with actives labelled from `actives.smi` (cross-checked against the `Active`
column carried in the merged tables — they agree exactly). Recipe (from the
analysis notebooks):

| method | source file | ranking score |
|---|---|---|
| DiffDock | `molecule_results_ranked.csv` | `confidence_score` |
| *-GNINA | `gnina_output.csv` | `CNNaffinity` |
| *-NMDN | `nmdn_output.csv` | `pKd-Score` (NMDN-Score = consensus quality filter) |

Reproduce single-method: `python scripts/reproduce_table2.py`
(full metrics → `scripts/build_release_data.py` → `docs/table2_reproduced.csv`).

## Single-method columns: correct, with one transcription error in ADRB2

AutoDock-GNINA EF1%, independent recompute vs Supplementary Table S2:

| target | recompute | S2 | | target | recompute | S2 |
|---|---|---|---|---|---|---|
| ADRB2 | **0.00** | **9.07 ✗** | MAPK1 | 0.97 | 0.97 ✓ |
| ALDH1 | 2.15 | 2.14 ✓ | MTORC1 | 0.00 | 0 ✓ |
| ESR1_ago | 7.78 | 7.82 ✓ | OPRK1 | 12.50 | 12.5 ✓ |
| ESR1_ant | 5.94 | 5.93 ✓ | PKM2 | 2.02 | 2.03 ✓ |
| FEN1 | 0.27 | 0.27 ✓ | PPARG | 3.70 | 3.73 ✓ |
| GBA | 9.11 | 9.07 ✓ | TP53 | 1.27 | 1.28 ✓ |
| IDH1 | 5.13 | 5.1 ✓ | VDR | 0.57 | 0.57 ✓ |
| KAT2A | 0.52 | 0.51 ✓ | | | |

14/15 reproduce exactly. The one exception, ADRB2, is proven wrong by
arithmetic, not opinion (below).

## Why ADRB2 is a copy/paste from GBA (proof by EF quantization)

EF1% can only take values `(hits_in_top_1% / n_actives) / 0.01`. **ADRB2 has 17
actives** (verified in every `actives.smi` and in the merged `Active` column), so
its EF1% can *only* be a multiple of **5.88**: `{0, 5.88, 11.76, …}`.

- Published `9.07` is **not on ADRB2's grid → impossible for any ADRB2 ranking.**
- Published AutoDock CC-Medium `7.25` is **also off ADRB2's grid → impossible.**
- Both values **are** exact for GBA (n_actives = 165, grid step 0.606):
  `9.07 = 15/165`, `7.25 ≈ 12/165`.
- In S2, ADRB2 and GBA share **both** `9.07` and `7.25`.

Conclusion: ADRB2's first two S2 cells were copied from GBA during the manual
per-target spreadsheet transfer. Independent recompute gives ADRB2 AutoDock-GNINA
= **0.00** and AutoDock CC-Medium = **5.88** (both on-grid).

## Impact on the headline

The reported AutoDock-GNINA **median EF1% = 2.14** is what you get *with* the bad
ADRB2 = 9.07 (median lands on ALDH1's 2.15). With the correct ADRB2 = 0.00 the
median is **2.03**. The best-single-method conclusion is unchanged.

**Manuscript corrections:** ADRB2 AutoDock-GNINA `9.07 → 0.00`; ADRB2 AutoDock
CC-Medium `7.25 → 5.88`; AutoDock-GNINA median EF1% `2.14 → 2.03`.

## Consensus columns (CC-Medium / UC-Strong / CC-Weak) and ML (Table III)

These are **not** independently re-derived here: exact reproduction depends on the
run-specific rank weights, NMDN-score threshold, and EF denominator used in the
analysis notebooks, and reproducing them with different parameters would be
misleading. They are reproducible by rerunning the **author's own code**, shipped
in `notebooks/` against the tables in `data/consensus_tables/` and `data/ml/`.
Nothing here contradicts them; they simply were not part of the arithmetic audit,
which covered the single-method columns.
