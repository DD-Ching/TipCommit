# Q2 Summary — TipCommit v0

**Date:** 2026-05-05
**Status:** Q2 v0 complete. Standing outputs frozen. Repo moves into
the second summary-mode pass.
**Predecessors:** [Q1 v0 summary](q1_summary.md) ·
[Q1.5b replication](q15b_replication.md) ·
[Q2 design](q2_design.md) · [Q2 MVP](q2_mvp.md) ·
[Q2.1 candidate-pool-only](q21_candidate_pool_only.md)

This note is the standing v0 summary of Q2 for a new researcher.
It pulls together the design-time question, the locked method, the
MVP's initial reading, the Q2.1 ambiguity-resolution pass, and the
final combined interpretation. **No new analysis.**

---

## 1. The Q2 question

> Within the He 2022 fetal lung epithelium, when cells leave the
> candidate distal-tip state and lean toward airway or alveolar
> fate, do they occupy a continuous gradient (smooth density between
> distal-tip and the two committed poles), a discrete switch
> (sparse gap between bipotent and committed states; bimodal
> commitment scores per stage), or a mixed pattern (gradient in one
> direction or one stage, discrete in another)?

Q2 was framed at design time as a **model-discrimination question**,
not a single-answer question. The MVP measured three diagnostics on
per-stage cell-state distributions and scored them against
pre-registered criteria for each pattern.

What Q2 v0 was **not** asking: when does commitment "happen" (= Q1,
already answered); what TF or pathway "drives" commitment (out of
scope, requires perturbation); how long an individual cell takes
to transit (out of scope, requires lineage tracing).

---

## 2. Substrate role split

### He 2022 — primary

Census `dataset_id == "3dc61ca1-ce40-46b6-8337-f27260fd9a03"` (the
He et al. 2022 fetal lung "All cells" release; same slice as Q1).
6 stages × ~10,300 epithelial cells; 10x 5′ v1 assay; 100% primary
data; AT1 / AT2 + airway sub-types separately annotated.

Why primary:
- **Full stage range (9–22 wpc).** Q2 needs both endpoints to
  discriminate gradient from switch. A substrate that misses the
  high-bipotent early window or the late committed window cannot
  test the question.
- **Reliable SOX2 / SOX9 detection** (Q1.5b showed Cao 2020's
  sci-RNA-seq3 detects SOX2 / SOX9 in ~10× fewer cells; Q2's
  bimodality test would be artefactually unimodal at low detection).
- **Both committed poles annotated separately.** Cao 2020 lacks
  AT1 / AT2 annotations.

### Cao 2020 — secondary, direction only

Census `dataset_id == "fa27492b-82ff-4ab7-ac61-0e2b184eee67"` (Cao
et al. 2020 "Survey of human embryonic development", 1M-cell
subset). 3 usable stages (15, 16, 17 wpc) × ~22,100 epithelial
cells.

Cao's role in Q2 was **strictly directional support** within its
tractable 15–17 wpc subset, not a quantitative replication. It did
not contribute to the headline pattern verdict, the gene
signatures, or the gradient-vs-switch criteria — per the Q2 design
and per Q1.5b's documented assay limits.

---

## 3. Locked multi-gene signatures

Three signatures, four genes each, anchored on canonical biology +
the markers that independently replicated in Q1.5b. **Locked at
design time; no signature tuning during MVP execution.**

| Signature | Genes | Anchor |
|---|---|---|
| `distal_tip_score` | SOX9, ID2, ETV5, TESC | Canonical SOX9+ tip; ETV5 replicated in 5/5 Cao stages |
| `airway_score` | SOX2, TP63, SCGB3A2, FOXJ1 | Airway TF + basal + secretory + multiciliated TF; SCGB3A2 replicated in 5/5 Cao stages |
| `alveolar_score` | SFTPC, SFTPB, AGER, SLC34A2 | AT2 surfactants + AT1 + AT2 transporter; SFTPC / SLC34A2 replicated in 5/5 Cao stages |

Notable design choices:
- **SOX2 in airway only** (despite tip-cell co-expression) so
  the distal-tip signature is orthogonal to airway for the
  bimodality test.
- **NKX2-1 excluded** from alveolar (pan-lung, would muddy
  airway-vs-alveolar contrast).
- **CLDN18 excluded** from alveolar (shared with stomach / GI).
- Scoring via `sc.tl.score_genes` with default control-set
  construction (`n_bins=25`, `ctrl_size=50`). No new dependencies.

The per-cell `commitment_score = airway_score − alveolar_score`
was the primary axis on which the bimodality test ran.

---

## 4. Q2 MVP — initial result (mixed at the all-epithelial level)

The MVP applied three diagnostics per stage in two views:
`all_epithelial` (the full ~10,300-cell set) and `exiting_tip`
(cells with `distal_tip_score` below the substrate-wide median —
a guard against attributing commitment shape to cells still
strongly tip-like). The user-required two-view refinement was the
only deviation from the design.

### 4a. He 2022 BC trajectory (both views agreed)

| Stage | BC all_epi | BC exiting_tip | Side of 0.555 |
|---|---:|---:|---|
| 9 wpc | 0.502 | 0.518 | gradient |
| 11 wpc | 0.470 | 0.415 | gradient |
| 15 wpc | 0.706 | 0.691 | **switch** |
| 18 wpc | 0.700 | 0.640 | **switch** |
| 20 wpc | 0.855 | 0.754 | **switch** |
| 22 wpc | 0.781 | 0.789 | **switch** |

Clear stage-dependent flip from gradient (BC < 0.555 at 9 + 11 wpc)
to switch (BC > 0.555 at 15 + 18 + 20 + 22 wpc), with the crossover
between 11 and 15 wpc — coincident with the candidate_bipotent
collapse identified in Q1 v0.

The joint-density grids confirmed the shape: 11 wpc was blob-like
(mass spread across the (airway, alveolar) plane); 18 + 22 wpc
showed two corners (committed alveolar and committed airway poles).

### 4b. Soft-switch caveat

The gap-zone fraction stayed high (47–90% in `all_epithelial`;
62–96% in `exiting_tip`) at every stage including the late, BC-
switch-like ones. A textbook hard switch would show gap-zone
< 10%; the observed values are 6–9× above that. The MVP labelled
this a **soft switch with persistent intermediate population**,
not a clean discrete step.

### 4c. Cao 2020 (secondary) reading

Cao BC was uniformly low (0.31–0.38) across 15–17 wpc — gradient-
like at every stage in both views, including the only directly
overlapping 15 wpc stage where He showed BC = 0.71. The
disagreement was attributed (per Q1.5b) to assay-level dropout
flattening Cao's score distributions, not to underlying biological
inconsistency. Cao supported gradient-like shape at 12–17 wpc but
was silent on the late-stage hardening that Q2 MVP suggested in
He 2022.

### 4d. Open ambiguity at MVP-time

Section 4a of [`notes/q2_mvp.md`](q2_mvp.md) flagged the central
unresolved question: was the late-stage hardening a **real
hardening of cell state**, or a **cell-composition artefact**
(more upstream-annotated AT1 / AT2 / basal / club / secretory /
multiciliated cells in the per-stage population at later stages
mechanically raising BC without anyone hardening within the
candidate pool)?

The exiting-tip view did not discriminate — the annotated
committed cells all have low `distal_tip_score` and survive the
exiting-tip filter; they are precisely the cells driving the
corner mass.

---

## 5. Q2.1 — candidate-pool-only resolution

The narrow Q2.1 follow-up restricted the same diagnostics to
`cell_type == "epithelial cell of lung"` only (the candidate pool
annotation; ~5,943 cells across all stages). Pure pandas
restructure of the Q2 per-cell scores; no Census re-fetch; no new
methods.

### 5a. Within-pool BC trajectory

| Stage | n in pool | BC all_in_pool | BC exit_in_pool | (Q2 MVP all_epi) |
|---|---:|---:|---:|---:|
| 9 wpc | 1,378 | **0.306** | 0.348 | 0.502 |
| 11 wpc | 1,409 | **0.332** | 0.363 | 0.470 |
| 15 wpc | 1,418 | **0.384** | 0.480 | 0.706 |
| 18 wpc | 1,254 | **0.333** | 0.353 | 0.700 |
| 20 wpc | 70 | 0.394 | 0.601 (n=40) | 0.855 |
| 22 wpc | 66 | 0.429 | 0.429 (n=50) | 0.781 |

**Every BC value with n ≥ 100 stays below the 0.555 bimodality
threshold in the all-in-pool view.** The 15 wpc within-pool BC
(0.38) is roughly half the all-epithelial value (0.71) — the
hardest direct test of where the late-stage switch was coming
from. The mixed pattern does **not** survive restriction to the
pool.

### 5b. Supporting in-pool signals

- **Gap-zone fraction stays ≥ 0.52 at every stage with n ≥ 100**
  in the all-in-pool view (0.82, 0.64, 0.52, 0.55 at 9, 11, 15,
  18 wpc). The candidate pool itself is dominated by intermediate-
  commitment cells throughout.
- **Co-commitment fraction peaks at 15 wpc inside the pool (0.34
  all view, 0.37 exit view).** High co-commitment is a textbook
  gradient signature (cells expressing both programs at once) —
  the opposite of what a switch would predict.

### 5c. Verdict

The Q2 MVP's late-stage switch was almost entirely driven by
inclusion of the upstream-annotated committed cells. The candidate
pool itself looks **gradient-like at every observable stage**.
The Q2 v0 ambiguity is resolved in favour of cell-composition
shift, not hardening.

---

## 6. Final combined interpretation

> Within the candidate pool, the bipotent → committed transition
> appears continuous (gradient-like) at every observable stage.
> The all-epithelial mixed pattern (gradient early, switch-like
> late) is the joint result of (1) a gradient-like in-pool
> distribution that stays gradient-like throughout development,
> plus (2) the Q1-documented depletion of that pool as cells
> transition into pre-existing committed-cell labels. Adding more
> committed-corner cells over time mechanically raises the all-
> cells bimodality coefficient without anyone hardening within
> the pool.

Restated in two sentences for a non-specialist:

- **The candidate bipotent pool keeps a smooth, intermediate-heavy
  airway-vs-alveolar distribution at every stage where we can
  measure it.**
- **What changes over time is how many cells are still in that
  pool, not the shape of the pool's commitment-score
  distribution.**

This is consistent with the textbook chain (NKX2-1 → SOX9+ tip →
bipotent → committed AT2 / AT1 or basal / secretory / ciliated)
without requiring any irreversible-switch event inside the
candidate pool itself.

The Cao 2020 secondary supports this reading: Cao's within-pool
BC at 15–17 wpc (0.31–0.38) is squarely in the gradient band,
matching He's within-pool BC at the same stages.

---

## 7. Explicit non-claims

The Q2 v0 result does **not** claim:

- **No lineage proof.** A snapshot expression-distribution test
  cannot prove that any individual cell did or did not transit
  through any state. Only lineage tracing can. `candidate_bipotent`
  remains a candidate label.
- **No causal proof.** Whether the airway-vs-alveolar choice is
  "driven by" SOX2, SOX9, FGF, WNT, BMP, or anything else is out
  of scope — Q2 only asks about the *shape* of the cell-state
  distribution, not its causes.
- **No regulatory proof.** Q2 does not test whether any TF binds
  any enhancer. The Gain-era public-data ChIP gap still applies.
- **No full trajectory inference / pseudotime.** Per design, no
  `sc.tl.dpt`, `sc.tl.diffmap`, scvelo, or RNA velocity. The
  signatures rank cells along committed scores at each stage; that
  is *not* pseudotime.
- **No claim that the signatures are exhaustive.** SOX9 + ID2 +
  ETV5 + TESC is one defensible distal-tip signature; alternative
  composites exist. v0 used a conservative locked set.
- **No claim of irreversibility at any week.** Even if some
  metric "hardens", v0 does not declare any cell irreversibly
  committed.
- **No claim about specific within-pool sub-structure.** Q2.1
  showed the pool stays gradient-like; v0 did not sub-cluster the
  pool.
- **No across-donor statistical inference.** Per-stage diagnostics
  are descriptive; with 10 donors across 6 stages, donor
  confounding is acknowledged.

---

## 8. Key limitations

1. **Late-stage candidate pool is too small for confident shape
   diagnostics.** The pool has 70 cells at 20 wpc and 66 at 22 wpc.
   Within-pool BC at those stages sits in the gradient band but
   with weaker statistical support than 9, 11, 15, 18 wpc.
2. **The candidate-pool annotation is itself a label.** "Epithelial
   cell of lung" was assigned by the upstream He 2022 annotation.
   If that annotation pre-classified some still-bipotent cells as
   committed, the within-pool reading would underestimate any
   real hardening. v0 trusts the upstream annotation, same as Q1.
3. **Two-gene definitions for transitioning sub-states (Q1.5b).**
   For backward compatibility with Q1, the SOX2-vs-SOX9 split
   inside the pool is still two-gene. Q2 used four-gene composites
   for the score-distribution test, which is more robust, but the
   transitioning_airway / transitioning_alveolar labels carried
   over from Q1 are not recomputed under the four-gene rule in v0.
4. **Cao 2020 cannot test the late-stage hardening question.** It
   stops at 17 wpc and its sci-RNA-seq3 sensitivity flattens the
   bimodality signal. The Q2.1 within-pool finding rests on He
   2022 alone.
5. **Single-pair signature comparison.** Only one airway signature
   and one alveolar signature were tested. A multi-signature
   ensemble (or alternative composites) might give different
   bimodality results. v0 is one defensible read, not the only
   one.
6. **No proximal-distal sub-axis.** The candidate pool was treated
   as one homogeneous group; if there is sub-structure
   (proximal-tip vs distal-tip), v0 cannot see it.
7. **Donor concentration per stage** carries over from Q1.

---

## 9. Where this leaves Q2

**The "shape" half of Q2 has a usable v0 answer:** within the
candidate pool, the bipotent-to-committed transition looks
continuous (gradient-like) at every observable stage. The
all-epithelial mixed pattern is a composition effect of pool
depletion, not within-pool hardening.

The complementary "what causes the choice" question — which is
distinct from gradient-vs-switch — remains entirely out of scope.
That is the kind of question that needs perturbation, lineage
tracing, or chromatin work; substrate-bound expression analysis
on its own cannot resolve it.

For the repo as a whole, see the updated [README](../README.md) for
the combined Q1 + Q1.5b + Q2 v0 standing result.
