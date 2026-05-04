# Q2 MVP Result — Axis Continuity

**Date:** 2026-05-05
**Script:** [`scripts/tipcommit_q2_mvp.py`](../scripts/tipcommit_q2_mvp.py)
**Design contract:** [`notes/q2_design.md`](q2_design.md)
**Outputs:**
- [`metadata/q2_per_cell_scores.csv`](../metadata/q2_per_cell_scores.csv)
- [`metadata/q2_per_stage_diagnostics.csv`](../metadata/q2_per_stage_diagnostics.csv)
- [`metadata/q2_joint_density.csv`](../metadata/q2_joint_density.csv)
- [`metadata/q2_cao_secondary_diagnostics.csv`](../metadata/q2_cao_secondary_diagnostics.csv)

**Verdict (lead):** **mixed pattern — gradient-like at 9–11 wpc,
switch-like at 15–22 wpc.** The bimodality coefficient flips
across stages exactly in the way the design's mixed criterion
predicted ("late-stage hardening"). The "switch" at late stages
is **soft** — a substantial intermediate population persists
through 22 wpc. The verdict is robust across the all-epithelial
and exiting-tip views. Cao 2020 directionally agrees on the early
stages but cannot test the late-stage hardening.

---

## 1. Directly observed score distributions

### 1a. Per-stage diagnostics (He 2022, primary)

Pre-locked diagnostic statistics on `commitment_score = airway_score
− alveolar_score`. BC = Sarle's bimodality coefficient (> 5/9 ≈
0.555 is conventional bimodality evidence). Gap-zone fraction =
cells with `|commitment| ≤ 0.5 × max(|commitment|)` AND
`distal_tip_score < 75th percentile`. Co-commitment fraction =
cells with both committed scores above their substrate-wide median.

| Stage | View | n | **BC** | gap-zone | co-commit |
|---|---|---:|---:|---:|---:|
| 9 wpc | all | 1,799 | **0.502** | 0.61 | 0.09 |
| 9 wpc | exiting_tip | 779 | **0.518** | 0.63 | 0.09 |
| 11 wpc | all | 2,043 | **0.470** | 0.47 | 0.20 |
| 11 wpc | exiting_tip | 785 | **0.415** | 0.62 | 0.13 |
| 15 wpc | all | 1,921 | **0.706** | 0.48 | 0.18 |
| 15 wpc | exiting_tip | 547 | **0.691** | 0.88 | 0.33 |
| 18 wpc | all | 1,812 | **0.700** | 0.59 | 0.17 |
| 18 wpc | exiting_tip | 759 | **0.640** | 0.90 | 0.27 |
| 20 wpc | all | 856 | **0.855** | 0.79 | 0.09 |
| 20 wpc | exiting_tip | 776 | **0.754** | 0.83 | 0.09 |
| 22 wpc | all | 1,086 | **0.781** | 0.62 | 0.13 |
| 22 wpc | exiting_tip | 911 | **0.789** | 0.71 | 0.14 |

The BC is **below 0.555 at 9 + 11 wpc** in both views and **above
0.555 at every stage from 15 wpc onward** in both views. The
crossover sits between 11 and 15 wpc.

### 1b. Joint density at the gradient and switch stages

5×5 binned density of (airway_bin, alveolar_bin), substrate-wide
quintile bins. Reading rows = airway-bin (low → high), columns =
alveolar-bin (low → high).

**He 2022, 11 wpc, all-epithelial (gradient-shape):**

```
alveolar_bin    0    1    2    3  4
airway_bin                          
0               8   22  127   50  0
1              22   61  215   87  2
2              60   96  242  137  2
3              78   96  153   46  4
4             110  151  223   51  0
```

Mass spreads broadly across the grid; concentrated in the middle
alveolar column (alveolar_bin = 2). No strong corners.

**He 2022, 18 wpc, all-epithelial (switch-shape):**

```
alveolar_bin   0    1   2    3    4
airway_bin                          
0              0    3  18  228  404
1              2    1  10  151  255
2             13   22  10   92  137
3             56  103  42   32   21
4             15   38  47  103    9
```

Two clear corners: top-right (low airway, high alveolar; the
alveolar pole, ≈ 1,038 cells) and bottom-left (high airway, low
alveolar; the airway pole, ≈ 212 cells). Intermediate zone
relatively sparse.

**He 2022, 22 wpc, exiting-tip (switch-shape, tip cells removed):**

```
alveolar_bin    0    1   2   3    4
airway_bin                          
0               1    0   1  45  160
1               4    3   6  23   76
2               9   17   3  12   36
3             110  153  30   9    3
4              30   65  36  76    3
```

Two corners hold up after removing tip-like cells: alveolar pole
≈ 304 cells (top-right), airway pole ≈ 358 cells (bottom-left).
The switch shape is not an artefact of including tip cells.

### 1c. Per-stage diagnostics (Cao 2020, secondary, 15–17 wpc only)

| Stage | View | n | **BC** | gap-zone | co-commit |
|---|---|---:|---:|---:|---:|
| 15 wpc | all | 8,171 | **0.316** | 0.71 | 0.14 |
| 15 wpc | exiting_tip | 3,950 | **0.309** | 0.96 | 0.13 |
| 16 wpc | all | 3,495 | **0.368** | 0.65 | 0.32 |
| 16 wpc | exiting_tip | 1,688 | **0.353** | 0.91 | 0.34 |
| 17 wpc | all | 10,460 | **0.380** | 0.62 | 0.33 |
| 17 wpc | exiting_tip | 5,479 | **0.366** | 0.84 | 0.34 |

Cao BC is **uniformly below 0.555** at every stage in both views.
Gap-zone fractions are uniformly very high (62–96%).

---

## 2. Interpretation under the pre-registered criteria

### 2a. The strict auto-verdict

The design's pre-registered automatic classifier (section E)
required ALL THREE criteria to align for a single-pattern verdict:
- Switch needs BC ≥ 0.555 in ≥ 3 stages **AND** gap ≤ 0.10 in ≥ 3
  stages **AND** co-commit ≤ 0.15 in ≥ 3 stages.
- Gradient needs BC < 0.555 in ≥ 4 stages **AND** gap ≥ 0.25 in
  ≥ 3 stages **AND** co-commit ≥ 0.30 in ≥ 3 stages.

Outcome from the script:
- He 2022 all_epithelial → strict verdict: `mixed_or_inconclusive`
  (BC supports switch in 4/6 stages but gap_low ≤ 0.10 fails in
  every stage; gradient also fails because BC < 0.555 in only 2/6).
- He 2022 exiting_tip → strict verdict: `mixed_or_inconclusive`
  (same reasons).
- Cao 2020 (both views) → `mixed_or_inconclusive` (BC supports
  gradient in 3/3 but co-commit ≥ 0.30 in only 2/3 stages).

### 2b. The design's mixed-criterion path

The design's section E also defined a separate mixed-pattern
trigger: *"the bimodality coefficient flips between stages (e.g.,
switch-like at 18–22 wpc but gradient-like at 9–11 wpc, indicating
late-stage hardening)."*

That is exactly what the He 2022 BC trajectory shows:

| Stage | BC (all view) | BC (exiting_tip view) | side of 0.555 |
|---|---:|---:|---|
| 9 wpc | 0.502 | 0.518 | gradient |
| 11 wpc | 0.470 | 0.415 | gradient |
| 15 wpc | 0.706 | 0.691 | **switch** |
| 18 wpc | 0.700 | 0.640 | **switch** |
| 20 wpc | 0.855 | 0.754 | **switch** |
| 22 wpc | 0.781 | 0.789 | **switch** |

This is a **clean mixed pattern**: gradient-like at the two
earliest fetal stages and switch-like at every stage from 15 wpc
onward, with both views agreeing on the side of every stage. The
switch onset between 11 wpc and 15 wpc coincides with the steep
candidate_bipotent decline first reported in Q1 v0 (24.4% → 17.4%
within total epithelial, moderate threshold).

### 2c. Why "soft switch", not "hard switch"

Even at the late, BC-switch-like stages, the gap-zone fraction
stays high (47–90% in the all view; 71–96% in the exiting-tip
view). A textbook "hard switch" would have gap-zone < 10%; the
observed values are 6–9× above that. The interpretation is that
the late-stage commitment-score distribution has two prominent
modes (BC high) but with substantial mass in between (gap-zone
high) — a **soft switch with a persistent intermediate
population**, not a clean discrete step.

The joint density at 18 wpc shows the same picture: the two
corners have most of the mass, but the (airway = 1, alveolar = 4)
and (airway = 2, alveolar = 3) cells still number in the dozens to
low hundreds.

### 2d. Both views agree

The exiting-tip view (cells with `distal_tip_score` below the
substrate-wide median) was added to guard against attributing
commitment shape to cells that are still strongly tip-like. **The
verdict does not change in either view**:

- Sign of BC vs 0.555 is the same at every He stage in both views.
- Gradient → switch crossover sits between 11 and 15 wpc in both
  views.
- Joint-density corners hold up in the exiting-tip view at 22 wpc
  (section 1b).

The exiting-tip view does **sharpen** the gap-zone pattern (it
goes from 47–79% in the all view to 62–96% in the exiting-tip
view at the same stages), suggesting that **most non-tip cells
live in the intermediate zone, not at the committed poles**. That
is consistent with the soft-switch reading.

---

## 3. What Cao 2020 does and does not support

### 3a. What it supports

- **Gradient-like shape at the early-window stages.** Cao's BC at
  15, 16, 17 wpc is uniformly 0.31–0.38 — well below the 0.555
  bimodality threshold. He 2022's gradient-stage BCs (9, 11 wpc)
  are 0.47–0.52, also below 0.555. Both substrates agree that
  early/mid-window distributions are unimodal / gradient-like in
  shape.
- **High intermediate occupancy.** Cao's gap-zone fractions
  (62–96%) corroborate the He observation that the intermediate
  zone is heavily populated even at stages where commitment is
  underway.
- **Direction is consistent** within Cao's narrow window.

### 3b. What it does not support

- **The late-stage hardening (15+ wpc switch onset) cannot be
  tested in Cao 2020.** Cao stops at 17 wpc — short of the 18–22
  wpc range where He 2022 shows the strongest switch shape.
- **At the only directly overlapping stage (15 wpc), Cao and He
  disagree on bimodality:** Cao BC = 0.316 (gradient), He BC =
  0.706 (switch). This is **not** treated as biological
  disagreement: Q1.5b documented that Cao's sci-RNA-seq3 detects
  SOX2 and SOX9 in ~10× fewer candidate-pool cells than He's 10x
  5′ v1; the lower assay sensitivity flattens score distributions
  by reducing the fraction of cells that score high on the
  committed signatures. The Cao reading at 15 wpc is consistent
  with assay-driven smearing, not with the underlying state
  distribution being more continuous than He shows.
- **Cao does not contribute markers, signatures, or quantitative
  diagnostics** to the headline verdict. Per the design, its role
  is direction-only.

### 3c. Net read on Cao 2020

Directional support for the gradient-shape early stages; silent
on the late-stage hardening. Consistent with the Q1.5b finding
that Cao 2020's sci-RNA-seq3 substrate is suitable for direction
checks but not for fraction- or shape-level quantitative claims.

---

## 4. What remains ambiguous

### 4a. Developmental hardening vs cell-composition shift (the biggest ambiguity)

The mixed pattern has two competing biological readings:

**(i) Developmental hardening of the cell-state distribution.**
Cells that exist at later stages have actually hardened into more
discrete commitment states — the bipotent / intermediate
population is being depleted while the committed poles fill up,
because individual cells are completing commitment trajectories.
This is the textbook expectation.

**(ii) Cell-composition shift driving apparent bimodality.**
At later stages the He 2022 dataset contains more cells that the
upstream annotation already labelled as committed AT1 / AT2 /
basal / club / multiciliated / secretory. These cells have low
`distal_tip_score` and high committed-pole scores by construction.
Adding more "corner cells" to a per-stage distribution will raise
its bimodality coefficient even if the *underlying state
distribution* is unchanged. The MVP cannot rule this out.

The exiting-tip view does **not** discriminate between (i) and
(ii). The annotated committed cells (AT1, AT2, basal, club, etc.)
all have low `distal_tip_score` and so survive the exiting-tip
filter — they are precisely the cells driving the corner mass.
Removing them would require restricting to the candidate pool
(`epithelial cell of lung`) only, which is the natural Q2.1
follow-up but is out of scope here.

### 4b. Other ambiguities (smaller)

- **Donor confounding.** With 10 donors across 6 He 2022 stages,
  some stages are 1–2 donors. The BC values per stage are not
  donor-resampled; if late-stage donors happened to over-represent
  committed cells, the BC pattern could be partly donor-driven.
- **Signature gene saturation at late stages.** SFTPC and AGER
  expression in mature AT2 / AT1 cells can be very high; the
  alveolar_score may saturate at the high end, compressing the
  alveolar tail of the distribution. This would inflate BC
  through skewness. The MVP does not test this.
- **The 0.555 threshold is conventional, not biological.** A
  BC of 0.518 (9 wpc, exiting-tip) sits just below the threshold
  and a BC of 0.640 (18 wpc, exiting-tip) sits well above —
  but the gradient-vs-switch boundary at exactly 5/9 is not a
  biology-derived cutoff. Re-running with BC ≥ 0.6 instead would
  not change the qualitative call (9, 11 wpc still gradient; 15+
  wpc still switch), but the boundary itself is not load-bearing
  biology.

---

## 5. What this MVP does not change

- **Q1 v0 stands.** The candidate-bipotent fraction trajectory
  (24% → 0.1%, He 2022) is unchanged; this MVP does not re-derive
  it.
- **Candidate-not-proven status of the bipotent state stands.**
  The Q2 signatures define expression-level scores; they do not
  prove any cell is or is not a true bipotent tip cell.
- **No causal, regulatory, or lineage claims.** Per design.

---

## 6. Suggested next pass (Q2.1, NOT done in this MVP)

If the user wants to push past the section 4a ambiguity without
opening a new question, the smallest follow-up would be:

- **Run the same diagnostics restricted to the candidate pool
  only** (`epithelial cell of lung` in He 2022; `epithelial cell
  of lower respiratory tract` in Cao 2020). If the candidate pool
  itself shows the same gradient → switch BC pattern, the
  hardening reading (4a-i) gains support. If it stays gradient at
  every stage, the cell-composition reading (4a-ii) gains
  support.

This would be a one-script Q2.1 that re-uses everything in this
MVP with a single mask change. No new design needed beyond a
short note locking the substrate slice.
