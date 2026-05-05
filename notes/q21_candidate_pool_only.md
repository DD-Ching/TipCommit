# Q2.1 — Candidate-Pool-Only Diagnostic Pass

**Date:** 2026-05-05
**Script:** [`scripts/tipcommit_q21_candidate_pool_only.py`](../scripts/tipcommit_q21_candidate_pool_only.py)
**Predecessor:** Q2 MVP — [`notes/q2_mvp.md`](q2_mvp.md)
**Output:** [`metadata/q21_candidate_pool_diagnostics.csv`](../metadata/q21_candidate_pool_diagnostics.csv)

This is a narrow ambiguity-test pass. **No new design, no Census
re-fetch, no new dependencies, no new biological scope.** It
re-runs the Q2 MVP's three diagnostics on the same per-cell scores
restricted to `cell_type == "epithelial cell of lung"` (He 2022's
candidate-pool annotation). Cao 2020 is not included — Q1.5b
already documented its assay-sensitivity unfitness for shape-level
claims.

---

## 1. The ambiguity being tested

From [`notes/q2_mvp.md`](q2_mvp.md) section 4a:

> The mixed pattern has two competing biological readings:
>
> (i) Developmental hardening of the cell-state distribution — cells
> at later stages have actually hardened into more discrete
> commitment states.
>
> (ii) Cell-composition shift driving apparent bimodality — at
> later stages the dataset contains more cells the upstream
> annotation already labelled as committed AT1 / AT2 / basal /
> club / multiciliated / secretory; adding these "corner cells"
> raises the bimodality coefficient mechanically without any
> underlying state-distribution change in the candidate pool
> itself.

Q2.1 directly tests the candidate-pool-only distribution: if (i),
the BC inside the candidate pool itself should show the
gradient → switch flip. If (ii), the candidate pool should look
gradient-like at every stage and the late-stage switch should
disappear.

---

## 2. Method (one paragraph)

Pure pandas restructure of `metadata/q2_per_cell_scores.csv`:
filter to He 2022 + `cell_type == "epithelial cell of lung"`,
recompute the three diagnostics (Sarle bimodality coefficient,
gap-zone fraction, co-commitment fraction) per stage in two views
(`all_in_pool`, `exiting_tip_in_pool`). Reference quantities
(distal-tip median for the exiting-tip cutoff, distal-tip 75th
percentile for the gap-zone cutoff, airway/alveolar medians for
co-commitment) are recomputed **within the candidate pool** —
that is the universe for this pass.

The diagnostic *definitions* (BC formula; gap-zone formula;
co-commitment formula) are unchanged. Only the input cell set
and the reference quantities change.

---

## 3. Headline — within-pool BC trajectory

| Stage | n in pool | n in exit-in-pool | **BC all_in_pool** | **BC exit_in_pool** |
|---|---:|---:|---:|---:|
| 9 wpc | 1,378 | 891 | **0.306** | **0.348** |
| 11 wpc | 1,409 | 554 | **0.332** | **0.363** |
| 15 wpc | 1,418 | 481 | **0.384** | **0.480** |
| 18 wpc | 1,254 | 664 | **0.333** | **0.353** |
| 20 wpc | 70 | 40 | **0.394** | **0.601** ⚠ small-n |
| 22 wpc | 66 | 50 | **0.429** | **0.429** ⚠ small-n |

**Every BC value in the all_in_pool view is below the bimodality
threshold (5/9 ≈ 0.555) at every stage.** In the exit_in_pool
view, only 20 wpc crosses the threshold (0.601) — and that is
based on **40 cells**, well into noise territory. Every stage with
a usable candidate-pool size (n ≥ 100) stays gradient-side at
both views.

---

## 4. Side-by-side BC comparison with Q2 MVP

| Stage | Q2 all_epi | Q2 exiting_tip | **Q2.1 all_in_pool** | **Q2.1 exit_in_pool** |
|---|---:|---:|---:|---:|
| 9 wpc | 0.502 | 0.518 | **0.306** | **0.348** |
| 11 wpc | 0.470 | 0.415 | **0.332** | **0.363** |
| 15 wpc | 0.706 ← switch | 0.691 ← switch | **0.384** | **0.480** |
| 18 wpc | 0.700 ← switch | 0.640 ← switch | **0.333** | **0.353** |
| 20 wpc | 0.855 ← switch | 0.754 ← switch | **0.394** | 0.601 (n=40) |
| 22 wpc | 0.781 ← switch | 0.789 ← switch | **0.429** | **0.429** |

(BC threshold for switch-like: 0.556)

The Q2 MVP showed BC crossing 0.555 between 11 and 15 wpc and
staying high through 22 wpc — the gradient → switch developmental
hardening reading. **That crossover does not survive the
restriction to the candidate pool.** Within-pool BC stays in the
0.30–0.48 band at every usable stage, well below the bimodality
threshold.

The 15 wpc within-pool BC (0.38) is the most direct comparison to
the Q2 MVP's first switch-like stage (0.71): a drop of nearly half.
The high BC at 15 wpc was being supplied almost entirely by the
already-committed-cell annotations.

---

## 5. Other diagnostics (full table)

Within-pool gap-zone fraction (cells in the intermediate
commitment zone with low distal-tip score) and co-commitment
fraction (cells above the 50th percentile on both committed
scores), per stage in both views:

| Stage | view | gap-zone | co-commit |
|---|---|---:|---:|
| 9 wpc | all_in_pool | 0.82 | 0.01 |
| 9 wpc | exit_in_pool | 0.94 | 0.00 |
| 11 wpc | all_in_pool | 0.64 | 0.13 |
| 11 wpc | exit_in_pool | 0.95 | 0.11 |
| 15 wpc | all_in_pool | 0.52 | 0.34 |
| 15 wpc | exit_in_pool | 0.82 | 0.37 |
| 18 wpc | all_in_pool | 0.55 | 0.25 |
| 18 wpc | exit_in_pool | 0.72 | 0.27 |
| 20 wpc | all_in_pool | 0.26 (n=70) | 0.14 |
| 20 wpc | exit_in_pool | 0.30 (n=40) | 0.13 |
| 22 wpc | all_in_pool | 0.18 (n=66) | 0.17 |
| 22 wpc | exit_in_pool | 0.24 (n=50) | 0.18 |

Two patterns worth flagging:

- **Gap-zone fraction stays high at every usable stage (≥ 0.52
  in the all view at 9, 11, 15, 18 wpc).** The "lots of
  intermediate cells" signal is not a Q2-MVP artefact; the
  candidate pool itself is dominated by intermediate-commitment
  cells throughout. The apparent late-stage gap-zone drop at 20–22
  wpc reflects the candidate pool's collapse to 70 / 66 cells, not
  a tightening of the distribution.
- **Co-commitment peaks at 15 wpc inside the pool (0.34 all-view,
  0.37 exit-view) and stays moderate through 18 wpc (0.25 / 0.27).**
  High co-commitment is a textbook gradient signature (cells
  expressing both programs at once). This is the opposite of what
  a switch reading would predict.

---

## 6. Verdict

**The mixed pattern does NOT survive restriction to the candidate
pool.** Within `epithelial cell of lung` cells alone, the
commitment-score distribution is **gradient-like at every stage
where n is large enough to read** (9, 11, 15, 18 wpc). BC stays
below the bimodality threshold; gap-zone stays high; co-commitment
peaks during the apparent transition window (15–18 wpc).

The Q2 MVP's "late-stage hardening" reading was real *as a
description of the all-epithelial set*, but it was driven by **cell-
composition shift** (the increasing fraction of upstream-annotated
committed cells at later stages), **not** by hardening of the
candidate pool itself. This is the section 4a-(ii) interpretation,
now strongly supported by direct evidence.

### What this changes about the Q2 reading

The honest combined reading from Q2 + Q2.1 is:

1. **Within the candidate pool, the bipotent → committed
   transition appears continuous (gradient-like) across the full
   9–18 wpc window.** Cells in the pool show smooth, intermediate-
   heavy commitment-score distributions; co-commitment is common.
2. **Cells leave the candidate pool to become annotated committed
   cells across the 11–22 wpc window** (Q1's main finding —
   24% → 0.1% within total epithelial). This depletion is the
   driver of the apparent population-level "switch."
3. **The all-epithelial mixed pattern is the joint result of
   (1) + (2):** a gradient-like in-pool distribution depleting
   into two pre-existing committed corners. Adding more committed-
   corner cells over time raises the all-cells BC mechanically
   without anyone hardening within the pool.

This is a meaningfully more conservative and more biologically
defensible read than "hardening." It is also consistent with the
Q1.5b Cao 2020 read (gradient-like across its 15–17 wpc window,
within the analog candidate pool).

---

## 7. Caveats

- **Late-stage candidate pool is too small for confident shape
  diagnostics.** At 20 wpc the pool has 70 cells; at 22 wpc, 66.
  The within-pool BC at those stages (0.39 / 0.43) sits in the
  gradient band but with weaker statistical support than the
  earlier stages. The exit-in-pool BC of 0.60 at 20 wpc is from
  40 cells — not load-bearing.
- **Candidate-pool annotation is itself a label.** "Epithelial
  cell of lung" was assigned by the upstream He 2022 annotation.
  If that annotation rule excluded some cells that are
  biologically still in the candidate pool but were
  pre-classified as committed, the within-pool reading would
  underestimate any hardening that actually happened within those
  cells. This is the same upstream-annotation-trust assumption
  Q1 v0 made.
- **The signatures themselves are unchanged.** Same four-gene
  composites for distal-tip, airway, alveolar. Q2's locked
  signature choices apply to Q2.1 verbatim; no signature tuning
  was performed.
- **No across-donor analysis.** Same donor-confounding caveat as
  Q1, Q1.5b, and Q2.

---

## 8. Net effect on Q2 readiness for summary

The main ambiguity flagged in [`notes/q2_mvp.md`](q2_mvp.md)
section 4a is **resolved in favour of cell-composition shift,
not hardening.** Q2 is now ready to be summarised cleanly with the
following corrected headline:

> Within the candidate pool, the bipotent → committed transition
> appears continuous (gradient-like) at every usable stage. The
> all-epithelial population shows a mixed gradient-then-switch
> pattern, but the late-stage "switch" is driven by cell-
> composition shift (more upstream-annotated committed cells over
> time) rather than by hardening within the candidate pool itself.

This is the same biological story as Q1 (the candidate pool
shrinks; cells move from candidate to committed labels) — Q2.1
shows that the cells *still* in the candidate pool at any given
stage continue to look gradient-like.

No further Q2 refinement passes are needed in this repo phase per
the user's instruction.
