# Q1 Key Findings — TipCommit v0

**Date:** 2026-05-05
**Full narrative:** [`notes/q1_summary.md`](q1_summary.md)

Concise enough for a new researcher to read in two minutes.

1. **The candidate SOX2 / SOX9 co-expressing fraction in human
   fetal lung epithelium declines strongly across development.**
   Within the He 2022 atlas's `epithelial cell of lung` candidate
   pool, the fraction at the moderate threshold (SOX2 and SOX9
   log-normalised expression both > 0.5) goes 21.8% → 24.4% → 17.4%
   → 4.7% → 0.2% → 0.1% across 9, 11, 15, 18, 20, 22 weeks
   post-conception. The 11 → 22 wpc decline is strictly monotonic
   and ~250×.

2. **The decline is threshold-robust.** Same shape at relaxed (>0),
   moderate (>0.5), and stringent (>1.0) cutoffs. Stricter
   thresholds compress the absolute fractions but do not change
   the trajectory.

3. **The decline is denominator-robust.** Same shape whether
   measured against total epithelial cells or against the candidate
   pool only (28.5% → 1.5% within-pool). The bipotent fraction
   collapses inside a pool that is itself shrinking — the trend
   is not an annotation artefact.

4. **There is a small early bump (9 → 11 wpc, +2.6 pp at moderate
   threshold) that should not be over-interpreted.** Two adjacent
   stages, small magnitude, donor-confounded. The post-11-wpc
   collapse is the load-bearing observation.

5. **Airway-shifted marker program is biologically coherent and
   stage-structured.** 15–18 wpc is dominated by secretory /
   apical / cytoskeletal genes (AGR2, AGR3, CAPS, ANXA1, EZR,
   PERP). 20–22 wpc is dominated by motile-cilia machinery
   (CFAP144/276, CIMAP1B, CIMIP1, SPMIP6, TSPAN1, DYNLL1).
   **SOX2 itself recurs as airway-shifted at 15 and 18 wpc** —
   internal validation that the upstream airway annotation is
   coherent.

6. **Alveolar-shifted marker program recovers canonical markers
   but is lower-magnitude.** AGER (AT1), SFTPB (AT2), CLDN18 (tight
   junction), MYL9, CLIC3, SMARCA5 recur across the 4 stages
   where DE could run. **NKX2-1** appears as alveolar-shifted at
   20 wpc. Ribosomal genes dominate the alveolar top-20s at every
   stage; the biological signal is genuine but smaller than the
   airway side's secretory / ciliated programs.

7. **Two stages have no testable airway-vs-alveolar DE.** 9 and
   11 wpc have **zero** annotated_alveolar cells in this dataset.
   The fraction trend is still reportable at those stages (no DE
   needed), but the marker comparison only covers 15, 18, 20, 22
   wpc.

8. **`candidate_bipotent` remains a candidate state, not a proven
   lineage.** SOX2 + SOX9 co-expression at scRNA-seq level is
   necessary but not sufficient for true tip-cell bipotency.
   v0 reports candidacy and a robust trend; functional /
   lineage-tracing proof is wet-lab and out of scope.
