# Findings — does multi-pose rescoring beat single-pose?

**Setup.** DiffDock pathway, 3 full-library LIT-PCBA targets, all 20 poses retained.
Single-pose (Arm A, published protocol) vs multi-pose rescoring of all 20 poses:
max-over-poses (B), mean-top-3 (C), NMDN pose-selector → pKd (D). Both rescorers
(GNINA, NMDN) run in their intended multi-pose mode — the thing Reviewer 3 asked for.
270,946 poses scored with NMDN + 258,000 with GNINA.

## Headline: EF1% unchanged in 5 of 6 target×scorer combos; one real GNINA gain

| Target | scorer | EF1% single (A) | EF1% multi-max (B) | EF1% mean-top3 (C) | verdict |
|--------|--------|-----------------|--------------------|--------------------|---------|
| TP53      | GNINA | 1.54  | 1.54  | 1.54  | no change |
| TP53      | NMDN  | 0.00  | 0.00  | 0.00  | no change |
| ESR1_ant  | GNINA | 2.22  | **3.33** | **4.44** | **improved (up to 2×)** |
| ESR1_ant  | NMDN  | 2.22  | 2.22  | 2.22  | no change |
| PPARG     | GNINA | 0.00  | 0.00  | 0.00  | no change |
| PPARG     | NMDN  | 20.83 | 20.83 | 20.83 | no change |

Two honest takeaways:

1. **NMDN — the method Reviewer 3 named — showed no EF1% gain on any of the three
   targets.** Multi-pose selection did not move early enrichment for NMDN at all
   (0→0, 2.22→2.22, 20.83→20.83).

2. **GNINA multi-pose helped on exactly one target (ESR1_ant), where mean-top-3
   doubled EF1% (2.22→4.44).** It did nothing on TP53 (1.54→1.54) and did not rescue
   PPARG (0→0). So the benefit is real but target-specific, not systematic.

## Broad-ranking metrics rise modestly under multi-pose

EF10%, ROC-AUC and BEDROC generally tick up with multi-pose (e.g. ESR1_ant NMDN
AUC 0.53→0.62; PPARG NMDN EF10% 3.75→4.58), i.e. extra poses refine mid-list order.
But except for ESR1_ant/GNINA this does not translate into the top 1% that decides a
screen.

## Arm D (NMDN pose-selector → pKd) is worse

Ranking by the pKd of the NMDN-selected pose collapses EF1% (0.00/1.11/0.00) — the
affinity readout is a poorer VS ranker than the NMDN score itself.

## Verdict / how to use in the rebuttal

We tested multi-pose rescoring directly, in the rescorers' intended mode, on three full
libraries. **Multi-pose did not systematically improve early enrichment: EF1% was
unchanged in 5 of 6 target×scorer cases, including all three NMDN targets.** The single
exception was GNINA on ESR1_ant, where mean-top-3 pose aggregation roughly doubled EF1%
(2.22→4.44). This supports two of the paper's claims at once: (i) single-pose is not a
crippling handicap for the top-1% comparison in general, and (ii) where rescoring gains
exist they are target- and method-specific rather than universal. We report the full
matrix and now state the single-pose limitation and this multi-pose control explicitly.

## Status
- NMDN arms (B/C/D): DONE, 3/3 targets.
- GNINA arms (B/C): DONE, 3/3 targets (fast parallel scorer, ~3.9 lig/s, GPU).
- All six target×scorer combinations complete. Raw numbers: `results/combined.csv`.

## Repro
```
scripts/assemble_poses.py <T>          # multi-model 20-pose SDFs (GNINA)
scripts/explode_poses.py  <T>          # per-pose SDFs (NMDN)
scripts/run_gnina_multipose.sh <T> 10  # GNINA, needs Docker (Linux CLI + /run/docker.sock)
scripts/nmdn_drive.py <T> --minutes 50 # NMDN, diffdock_nmdn env, resumable
scripts/aggregate_gnina.py <T> ; scripts/aggregate_nmdn.py <T> ; scripts/combine.py
```
