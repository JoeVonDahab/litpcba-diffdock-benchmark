# Draft reviewer response — multi-pose experiment

Covers Reviewer 3 points 1–2 and the related single-pose concerns from Reviewers 1 and 2.
Numbers from the multi-pose experiment (DiffDock pathway; 3 full-library LIT-PCBA targets
TP53/ESR1_ant/PPARG; all 20 poses rescored by GNINA and NMDN). Full matrix in
`results/combined.csv` and Supplementary Table SX.

---

**Reviewer 3, comments 1–2 (single-pose biases against rescorers, esp. DiffDock-NMDN,
which is designed to exploit multiple poses).**

We thank the reviewer and agree that single-pose evaluation does not exercise rescorers
in their intended operating mode. To measure the practical impact, we re-ran the DiffDock
pathway in full multi-pose mode on three full-library LIT-PCBA targets screened without
inactive subsampling (TP53, ESR1_ant, PPARG; 3,784–5,093 ligands). For every ligand we
rescored all 20 DiffDock poses with both GNINA and NMDN and selected the best-rescored
pose (maximum score across poses, and the mean of the top three), i.e. the standard
multi-pose protocol and the mode for which DiffDock-NMDN was designed.

Multi-pose rescoring did not systematically improve early enrichment. EF1% was unchanged
versus single-pose in five of the six target×scorer combinations, including all three
NMDN targets (NMDN EF1% 0.00→0.00, 2.22→2.22, 20.83→20.83 for TP53/ESR1_ant/PPARG; GNINA
1.54→1.54 on TP53 and 0.00→0.00 on PPARG). The single exception was GNINA on ESR1_ant,
where aggregating over poses increased EF1% from 2.22 to 3.33 (best pose) and 4.44
(mean of top three). Broader-ranking metrics (EF10%, ROC-AUC, BEDROC) rose modestly under
multi-pose across targets, indicating that additional poses refine mid-list ordering, but
this generally did not propagate to the top 1% that governs a prospective screen.

We have added the full matrix as Supplementary Table SX and revised the Methods and
Discussion to (i) state explicitly that single-pose evaluation does not represent the
intended multi-pose mode of rescorers such as DiffDock-NMDN, (ii) report that, tested in
that mode on full libraries, multi-pose rescoring improved early enrichment in only one
of six cases (GNINA/ESR1_ant), and (iii) moderate our conclusions accordingly: where
rescoring gains occur they are target- and method-specific rather than systematic, which
is consistent with the paper's central finding. This also aligns with the original
DiffDock-NMDN report, whose gains derived from far larger pose ensembles together with
target-specific tuning, neither of which applies in a fixed-protocol unbiased screen at
this scale.

---

**One-paragraph version (if space-limited):**

To address the concern that single-pose evaluation handicaps rescorers, we re-ran the
DiffDock pathway in full multi-pose mode (all 20 poses rescored by GNINA and NMDN,
best-rescored pose retained) on three full-library LIT-PCBA targets. Early enrichment did
not improve systematically: EF1% was unchanged in five of six target×scorer cases,
including all three NMDN targets, with the single exception of GNINA on ESR1_ant (EF1%
2.22→4.44 with mean-of-top-3 pose aggregation). We report the full matrix as
Supplementary Table SX and have explicitly flagged the single-pose limitation and this
multi-pose control in the revised Methods and Discussion.
