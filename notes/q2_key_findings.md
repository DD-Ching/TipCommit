# Q2 Key Findings — TipCommit v0

**Date:** 2026-05-05
**Full narrative:** [`notes/q2_summary.md`](q2_summary.md)

Two-minute read for a new researcher landing on the Q2 part of
the repo cold.

1. **Q2 is a model-discrimination question, not a single-answer
   question.** Pre-registered three patterns (gradient / switch /
   mixed) and three diagnostics (Sarle bimodality coefficient,
   gap-zone fraction, co-commitment fraction) on a per-cell
   `commitment_score = airway_score − alveolar_score`, scored on
   He 2022 fetal lung (primary) and Cao 2020 (secondary, direction
   only).

2. **Three multi-gene signatures (locked at design time, four
   genes each):** distal_tip = SOX9 + ID2 + ETV5 + TESC; airway =
   SOX2 + TP63 + SCGB3A2 + FOXJ1; alveolar = SFTPC + SFTPB + AGER
   + SLC34A2. Anchored on canonical biology + the markers that
   independently replicated in Q1.5b. Scored via `sc.tl.score_genes`.

3. **At the all-epithelial level the He 2022 bimodality
   coefficient flips from gradient-like at 9–11 wpc (BC ≈ 0.50)
   to switch-like at 15–22 wpc (BC ≈ 0.70–0.86)**, crossover
   between 11 and 15 wpc — coincident with Q1's candidate-bipotent
   collapse window. The exiting-tip view (cells below the
   distal-tip-score median) agrees on stage-by-stage shape.

4. **The "switch" at late stages is soft, not hard.** Even in the
   BC-switch-like stages, the gap-zone fraction (intermediate-
   commitment cells with low tip-character) stays at 47–90% — a
   textbook hard switch would show < 10%. A substantial
   intermediate population persists through 22 wpc.

5. **Q2.1 (the candidate-pool-only follow-up) overturned the
   late-stage switch reading.** Restricted to `epithelial cell of
   lung` cells alone, the within-pool BC stays at 0.31–0.43 at
   every stage with n ≥ 100 — uniformly **below** the 0.555
   bimodality threshold. The 15 wpc within-pool BC (0.38) is
   roughly half the all-epithelial value (0.71).

6. **The all-epithelial mixed pattern was driven by cell-
   composition shift, not by hardening within the candidate
   pool.** As cells transition out of the candidate pool into the
   pre-existing AT1 / AT2 / basal / club / multiciliated /
   secretory annotations (Q1's main finding: 24% → 0.1%), adding
   more committed-corner cells over time mechanically raises the
   all-cells BC. The candidate pool itself stays gradient-like
   throughout — including with high co-commitment (cells expressing
   both programs at once, peaking at 34% within-pool at 15 wpc).

7. **Cao 2020 directionally supports the gradient-like in-pool
   reading at 15–17 wpc** (within-pool BC analog 0.31–0.38) but is
   silent on the late-stage stages He 2022 covers (Cao stops at
   17 wpc, and Q1.5b documented its sci-RNA-seq3 sensitivity gap).

8. **`candidate_bipotent` remains a candidate state, not a
   proven lineage.** The shape-of-distribution test does not
   prove any cell is or is not a true bipotent tip cell.
