# Q2 Design — Axis Continuity (gradient vs switch)

**Date:** 2026-05-05
**Status:** design contract for Q2. **No implementation in this turn.**
**Predecessors:**
- Q1 v0 — [`notes/q1_summary.md`](q1_summary.md)
- Q1.1 refinement — [`notes/q1_refinement.md`](q1_refinement.md)
- Q1.5b replication — [`notes/q15b_replication.md`](q15b_replication.md)

This is a **single-question design**: the bipotent → committed
transition along the proximal-distal / airway-alveolar axis — is
the cell-state distribution continuous, discrete, or mixed?
Everything else is out of scope for Q2 v0.

---

## A. The exact Q2 question

> Within the He 2022 fetal lung epithelial atlas, when cells leave
> the candidate distal-tip state and lean toward airway or alveolar
> fate, do they occupy a **continuous gradient** (cell-state space
> filled with a smooth density between distal-tip and the two
> committed poles), a **discrete switch** (sparsely-populated gap
> between the bipotent state and the committed states, with a
> bimodal distribution of commitment scores per stage), or a
> **mixed pattern** (gradient in one direction or one stage,
> discrete in another)?

Q2 is **a model-discrimination question, not a single-answer
question.** The MVP measures three diagnostics on the cell-state
distributions and scores them against pre-registered criteria for
each of the three patterns.

What Q2 is **not** asking:
- "When does commitment happen?" — that is Q1, already answered.
- "What gene is the master regulator?" — out of scope; needs
  perturbation.
- "How long does a transitioning cell take to commit?" — needs
  lineage tracing or live imaging; out of scope.

---

## B. Primary substrate (He 2022) and secondary substrate (Cao 2020)

### He 2022 — primary

Same Census slice as Q1 v0:

```python
value_filter = (
    "tissue_general == 'lung' and "
    "dataset_id == '3dc61ca1-ce40-46b6-8337-f27260fd9a03'"
)
```

Reasons He 2022 is the right primary substrate for Q2:

1. **Full stage range (9–22 wpc).** The bipotent → committed
   transition window spans 11 → 18 wpc in v0. Q2 needs both
   endpoints to discriminate gradient from switch — a substrate
   that misses the early high-bipotent stages or the late
   committed stages cannot test the question.
2. **Reliable SOX2 / SOX9 + tip-marker detection** (10x 5′ v1).
   Q1.5b showed Cao 2020's sci-RNA-seq3 detects SOX2 / SOX9 in
   ~10× fewer cells than He 2022. Q2's diagnostics (bimodality,
   gap density) are sensitive to detection rate; running them on
   a low-sensitivity assay would give artefactually unimodal
   distributions for non-biological reasons.
3. **Annotated AT1 / AT2 + airway sub-types.** Q2's gene
   signatures are anchored on canonical committed-cell biology;
   He 2022 is the only available substrate where both committed
   poles are separately annotated (Cao 2020 lacks AT1 / AT2).
4. **Single assay, 100% primary data.** No cross-assay batch
   correction is needed inside the substrate.

### Cao 2020 — secondary, directional support only

Same Census slice as Q1.5b:

```python
value_filter = (
    "tissue_general == 'lung' and "
    "dataset_id == 'fa27492b-82ff-4ab7-ac61-0e2b184eee67'"
)
```

Cao 2020's role in Q2 is **strictly directional support** within
its tractable 15–17 wpc subset, not a quantitative replication.
Specifically: if the He 2022 bimodality / gap diagnostics point to
a clear pattern (gradient vs switch), do the same diagnostics on
Cao 2020 cells (computed at 15, 16, and 17 wpc only — the stages
where Cao has > 3,000 candidate-pool cells) point in the **same
direction**?

Cao 2020 is **not** used to:
- compute headline pattern verdicts,
- compute the discrete-zone fraction or gap-density quantitatively,
- compute a per-stage trajectory (its stage range is too narrow),
- contribute markers to the gene signatures (the signatures are
  pre-locked from prior literature + Q1 + Q1.5b).

If Cao 2020's diagnostics disagree with He 2022 on direction, that
disagreement is reported and attributed to assay sensitivity (per
Q1.5b's documented limitation), not used to reverse the He 2022
verdict.

---

## C. State signatures (multi-gene, conservative, pre-locked)

Three signatures, four genes each. Anchored on canonical biology
plus the markers that **independently replicated** in Q1.5b.
**Locked here in this design.** No tuning during MVP execution.

### C.1. Distal-tip signature (`distal_tip_score`)

Genes co-expressed in the canonical SOX9+ distal-tip / bipotent
progenitor in human fetal lung.

| Gene | Why included |
|---|---|
| **SOX9** | Canonical distal-tip TF (textbook). Defines the tip in mouse and human. |
| **ID2** | Tip / proliferating-progenitor TF in fetal lung; established marker in human and mouse fetal-lung literature. |
| **ETV5** | FGFR2-downstream TF. **Replicated cleanly at top ranks in 5/5 Cao 2020 stages** (Q1.5b) — strongest non-trivial recovery in the alveolar/distal direction. |
| **TESC** | Tescalcin; distal-tip marker described in human fetal lung scRNA-seq (Nikolic et al. 2017 era). |

Notes on what is **deliberately excluded**:
- SOX2 is excluded from this signature (despite tip cells
  co-expressing it) so that the distal-tip signature is
  orthogonal to the airway signature for the bimodality test.
- ID2's expression is partially shared with proliferating cells
  generally — interpreted alongside cell-cycle scoring is a v1+
  refinement; v0 takes ID2's signal at face value.

### C.2. Airway-leaning signature (`airway_score`)

Genes marking the proximal / airway program. Mixes one TF, one
secretory marker, one basal marker, one ciliated TF.

| Gene | Why included |
|---|---|
| **SOX2** | Canonical airway TF (textbook); appeared as airway-shifted in 2 of 4 Q1 stages and is the upstream regulator. |
| **TP63** | Basal-cell TF; canonical proximal/airway lineage marker. |
| **SCGB3A2** | Secretory / club-lineage marker. **Replicated at rank 2 in 5/5 Cao 2020 stages** (Q1.5b) — strongest non-trivial recovery in the airway direction. Known SOX2 target. |
| **FOXJ1** | Multiciliated-cell TF. The 20–22 wpc Q1 markers were dominated by motile-cilia genes; FOXJ1 is the upstream TF for that program. |

### C.3. Alveolar-leaning signature (`alveolar_score`)

Genes marking the distal / alveolar program. Mixes the two AT2
canonical genes that survived independent replication, the AT1
canonical, and an alveolar tight-junction marker.

| Gene | Why included |
|---|---|
| **SFTPC** | Canonical AT2 surfactant protein. **Replicated at rank 2 in 5/5 Cao 2020 stages** (Q1.5b). |
| **SFTPB** | Canonical AT2 / alveolar surfactant. Recurrent at 2/4 stages in He 2022 (Q1). |
| **AGER** | Canonical AT1 marker (RAGE). Recurrent at 4/4 stages in He 2022 (Q1). |
| **SLC34A2** | AT2 phosphate transporter (NaPi-IIb). **Replicated at top ranks in 5/5 Cao 2020 stages** (Q1.5b). |

NKX2-1 is **not** in this signature: it is pan-lung-epithelial
(expressed in both alveolar and airway lineages from endoderm
specification), so its inclusion would muddy the alveolar-vs-airway
contrast. Recurrence of NKX2-1 at 20 wpc as alveolar-shifted in Q1
is acknowledged but not used as a defining gene here.

CLDN18 is **not** in this signature despite Q1's recurrence at
2/4 stages: it is shared with stomach / GI epithelium and might
muddy the alveolar specificity at the broader-scope statistical
test. The four genes above are sharper for alveolar fate.

### Scoring method

`scanpy.tl.score_genes` with default control-set construction
(`n_bins=25`, `ctrl_size=50`), one call per signature. The output
is a continuous per-cell score (mean expression of signature genes
minus mean expression of a matched random control set). No
binarisation at the cell level. **No new packages.**

---

## D. The smallest useful Q2 MVP

A single new script `scripts/tipcommit_q2_mvp.py` that:

1. Re-fetches the He 2022 epithelial slice exactly as Q1 v0 did.
2. Re-uses Q1's normalisation pipeline.
3. Computes **three per-cell signature scores** via
   `sc.tl.score_genes` using the locked gene sets in section C.
4. Computes a **commitment score** per cell:
   `commitment = airway_score − alveolar_score`. Positive =
   airway-leaning; negative = alveolar-leaning; near-zero =
   uncommitted (or doubly-committed — distinguished by
   `distal_tip_score`).
5. Per stage, computes three diagnostic statistics on the
   `commitment` distribution:
   - **bimodality coefficient** (BC) — Sarle's BC =
     `(skew² + 1) / (kurt + 3(n-1)²/((n-2)(n-3)))`. BC > 5/9
     (≈ 0.555) is conventional evidence of bimodality. Computable
     with `scipy.stats` (already a transitive dep).
   - **gap-zone fraction** — fraction of cells whose
     `|commitment| ≤ 0.5 × max(|commitment|)` AND
     `distal_tip_score` is below its 75th percentile. Cells in
     the "intermediate but not distal" zone.
   - **co-commitment fraction** — fraction of cells with both
     `airway_score` and `alveolar_score` above their respective
     50th percentiles. High co-commitment = gradient (cells share
     both programs); low co-commitment = switch (mutually
     exclusive).
6. Per stage, also reports the joint 2D density (airway_score,
   alveolar_score) summarised as a 5×5 binned count matrix —
   small enough to inspect in a CSV without plotting.
7. Re-runs the same per-cell scoring on Cao 2020 (the Q1.5b
   substrate) at 15–17 wpc only, and reports the same three
   diagnostics for direction comparison only.
8. Writes outputs:

| Output | Shape | Purpose |
|---|---|---|
| `metadata/q2_per_cell_scores.csv` | one row per epithelial cell × 3 scores + cell_type + stage + donor | full data for downstream re-analysis |
| `metadata/q2_per_stage_diagnostics.csv` | one row per (stage, substrate) | bimodality coefficient, gap-zone fraction, co-commitment fraction, n_cells |
| `metadata/q2_joint_density.csv` | one row per (stage, airway_bin, alveolar_bin) × count | 2D density for visual inspection |
| `metadata/q2_cao_secondary_diagnostics.csv` | one row per Cao stage | direction-only check |

No plotting. No clustering. No pseudotime. No multi-batch
integration. Single script, ~300 LoC ceiling.

---

## E. Gradient vs switch criteria — pre-registered

Decision rules locked here before any data is seen. **The MVP does
not get to tune these after observing results.**

Let *He stages used* = 5 stages where the candidate pool has
≥ 100 cells (Q1 v0 found this is 9, 11, 15, 18 wpc — 20 and 22 wpc
are too small at moderate threshold for reliable bimodality
estimates and will be reported but not counted toward the verdict).
Replace if Q2's actual usable-stage set differs after the run.

### Switch-like verdict requires ALL THREE:

1. **Bimodality coefficient ≥ 0.555 in ≥ 3 of the usable stages.**
2. **Gap-zone fraction ≤ 10% in ≥ 3 of the usable stages.** (Few
   cells in the in-between zone.)
3. **Co-commitment fraction ≤ 15% in ≥ 3 of the usable stages.**
   (Cells largely commit to one program at a time.)

### Gradient-like verdict requires ALL THREE:

1. **Bimodality coefficient < 0.555 in ≥ 4 of the usable stages.**
2. **Gap-zone fraction ≥ 25% in ≥ 3 of the usable stages.** (Many
   cells in the in-between zone.)
3. **Co-commitment fraction ≥ 30% in ≥ 3 of the usable stages.**
   (Cells routinely score on both programs simultaneously.)

### Mixed verdict if EITHER:

- The bimodality coefficient flips between stages (e.g., switch-
  like at 18–22 wpc but gradient-like at 9–11 wpc, indicating
  late-stage hardening), AND the other two diagnostics agree at
  the relevant stages.
- One signature direction is bimodal but the other is unimodal at
  the same stage (e.g., airway commitment is switch-like; alveolar
  is gradient).

### Inconclusive verdict if:

- None of the above patterns hold cleanly. Reported as honest
  inconclusive; Q2 would then need a v1+ refinement (multi-gene
  signature expansion, alternative scoring method, or a different
  substrate).

### Cao 2020 supports / disagrees

Cao 2020 supports the He 2022 verdict if Cao's bimodality
coefficient at 15, 16, 17 wpc trends in the same direction
(higher / lower) as He 2022 at the comparable stages. If Cao
disagrees, the He verdict still stands but is flagged as
substrate-dependent; the Cao disagreement is documented and
attributed to assay sensitivity per Q1.5b.

---

## F. Explicit non-claims

The Q2 MVP will NOT claim:

- **No lineage proof.** A snapshot expression-distribution test
  cannot prove that any individual cell did or did not transit
  through any state. Only lineage tracing can.
- **No causal proof.** Whether the airway-vs-alveolar choice is
  "driven by" SOX2, SOX9, FGF, WNT, BMP, or anything else is out
  of scope — Q2 only asks about the *shape* of the cell-state
  distribution, not its causes.
- **No regulatory proof.** Q2 does not test whether any TF binds
  any enhancer. The chromatin-substrate gap from Gain still
  applies; Q2 makes no claim about TF–target relationships.
- **No full trajectory inference.** The signatures rank cells
  along committed scores at each stage; this is **not pseudotime**.
  No `sc.tl.dpt`, no `sc.tl.diffmap`, no `scvelo`. If a v1+ Q2
  follow-up wants pseudotime, it gets its own design note.
- **No claim that the signatures are exhaustive.** SOX9 + ID2 +
  ETV5 + TESC is one defensible distal-tip signature; alternative
  composites exist. The signature set is conservative and locked
  for v0; expansion is a v1 move.
- **No claim about specific timing of irreversibility.** Even if
  the verdict is "switch-like by 18 wpc", the MVP does not declare
  any individual cell irreversibly committed at any week.
- **No claim about within-pool sub-structure.** Q1.5a (sub-
  clustering the candidate pool) was deferred; Q2 does not run
  sub-clustering. The signature scores are computed per cell
  independently.
- **No across-donor inference.** With 10 donors across 6 stages
  (He) and 11 donors across 5 stages (Cao), some stages remain
  donor-confounded. Per-stage diagnostics are descriptive.
- **No claim that "gradient" rules out future commitment.**
  Observing a smooth distribution at a given snapshot does not
  imply cells stay smoothly distributed forever; only that the
  snapshot at that stage is smooth.

---

## G. Reuse of existing TipCommit logic

The Q2 MVP reuses, with no edits:

- The `fetch_anndata`, `filter_to_epithelial`, `normalise`,
  `get_gene_expression` helpers from
  `scripts/tipcommit_q1_mvp.py` (copy-paste; no shared module,
  per the project's "no generic framework" rule).
- The same Census version pin (`2025-11-08`).
- The same epithelial cell-type filter from Q1 v0.
- The Cao 2020 substrate constants from
  `scripts/tipcommit_q15b_replication.py`.

Q2 adds, and only adds:

- The three locked gene signatures (section C).
- One `sc.tl.score_genes` call per signature.
- Three per-stage diagnostic functions (BC, gap-zone fraction,
  co-commitment fraction).
- One 2D-histogram aggregator for the joint density CSV.

No new external packages.

---

## H. Out of scope for Q2 v0

- Pseudotime, RNA velocity, scvelo, scvi-tools.
- Sub-clustering the candidate pool (= Q1.5a).
- Cross-dataset integration / batch correction.
- Per-donor variance decomposition.
- Multi-gene signature expansion beyond the 4-gene-per-state
  composites.
- Functional / regulatory / lineage validation.
- Spatial transcriptomics cross-check.
- Plotting infrastructure (the joint-density CSV is the figure-
  ready substrate).

---

## I. Implementation prerequisites (for the next turn)

Before writing `scripts/tipcommit_q2_mvp.py`, confirm:

1. All 12 signature genes (SOX9, ID2, ETV5, TESC, SOX2, TP63,
   SCGB3A2, FOXJ1, SFTPC, SFTPB, AGER, SLC34A2) are present in
   the Census `var` annotations with `feature_name` matches. Q1.5b
   already showed SCGB3A2, SFTPC, ETV5, SLC34A2, SOX2, SOX9 are
   present; the remaining six (ID2, TESC, TP63, FOXJ1, SFTPB,
   AGER) are very common gene symbols and almost certainly
   present, but the MVP script should fail loudly if any are
   missing rather than silently dropping them.
2. `scipy.stats.skew` and `scipy.stats.kurtosis` (already a
   transitive dep) work on the score arrays as expected.
3. `sc.tl.score_genes` runs on the He 2022 AnnData without
   requiring extra preprocessing beyond the Q1 normalisation.

If any prerequisite fails, the design is revisited before
implementation. The expected case is all three pass.
