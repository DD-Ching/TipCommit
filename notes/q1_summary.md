# Q1 Summary — TipCommit v0

**Date:** 2026-05-05
**Status:** Q1 v0 result complete. Standing outputs frozen. Repo
moves into summary phase.
**Predecessor:** [Gain v0.1.0](https://github.com/DD-Ching/Gain/releases/tag/v0.1.0)
(evidence-audit / substrate-gap repo; the chain stopped at the
public-data ChIP gap, and TipCommit picked up the expression-
correlation route).

This note is the standing v0 summary for a new researcher. It pulls
together what was actually observed, what is robust, what stays
explicitly hypothetical, and where the analysis ran out of substrate.
Nothing here is new analysis — it summarises the MVP
([`scripts/tipcommit_q1_mvp.py`](../scripts/tipcommit_q1_mvp.py)) and
the Q1.1 refinement
([`scripts/tipcommit_q1_refinement.py`](../scripts/tipcommit_q1_refinement.py)).

---

## 1. The Q1 question

> When and how does the bipotent SOX2 / SOX9 co-expressing tip cell
> commit to airway vs alveolar fate in **human** fetal lung?

Verbatim from [`notes/evidence_map.md`](evidence_map.md).
TipCommit v0 only addresses the "when" half — the per-stage shape
of the candidate co-expressing fraction. The "how" half (continuous
gradient vs discrete switch — Q2) is explicitly deferred.

## 2. Dataset

**CELLxGENE Census** ("stable" LTS), human fetal lung, dataset
`3dc61ca1-ce40-46b6-8337-f27260fd9a03` — the integrated "All cells"
release of the He et al. 2022 collection, *"A human fetal lung cell
atlas uncovers proximal-distal gradients of differentiation and
key regulators of epithelial fates."*

| Attribute | Value |
|---|---|
| Cells (fetal-stage subset) | 62,759 |
| Of which epithelial (9 cell types) | 10,297 |
| Donors | 10 |
| Assay | 10x 5′ v1 (single platform) |
| Primary data fraction | 100% |
| Developmental stages | 9, 11, 15, 18, 20, 22 wpc |

The single-assay / single-collection design means no cross-platform
batch correction was needed inside v0. Substrate inventory is in
[`notes/he2022_census_inventory.md`](he2022_census_inventory.md).

## 3. Cell-state definitions

### 3a. Candidate pool

The cell-type annotation **`epithelial cell of lung`** (5,943 cells
in this run; ~5,595 in the inventory's fetal-only count) is treated
as the **candidate starting pool**. It is the residual epithelial
annotation in the upstream paper — cells the original clustering
did not place into a committed (basal / club / multiciliated /
secretory / AT1 / AT2 / squamous) bin.

This pool is **a candidate substrate, not a proven bipotent
population**. The MVP's role is to ask whether *within* this pool,
SOX2/SOX9 co-expression resolves into a sensible developmental
trajectory.

### 3b. Strict committed states (trust the upstream annotation)

| State | Cell types included |
|---|---|
| `annotated_airway` | basal, club, lung secretory, lung multiciliated, respiratory tract epithelial, squamous epithelial |
| `annotated_alveolar` | pulmonary alveolar type 1 (AT1), pulmonary alveolar type 2 (AT2) |

These labels are taken directly from the He 2022 annotation; the
MVP does not re-cluster or re-annotate.

### 3c. Within-pool sub-classification (expression-based)

After Scanpy normalisation (`sc.pp.normalize_total(target_sum=1e4)`
+ `sc.pp.log1p`), each candidate-pool cell is assigned to one of:

| Sub-state | SOX2 (lognorm) | SOX9 (lognorm) |
|---|---|---|
| `candidate_bipotent` | > threshold | > threshold |
| `transitioning_airway` | > threshold | ≤ threshold |
| `transitioning_alveolar` | ≤ threshold | > threshold |
| `undefined` | ≤ threshold | ≤ threshold |

Threshold sweep: relaxed > 0, **moderate > 0.5 (default headline)**,
stringent > 1.0. The moderate threshold was chosen after the user
flagged that the >0 threshold was too permissive for biological
interpretation.

## 4. The per-stage trend (headline result)

`candidate_bipotent` cells as a fraction of total epithelial cells,
moderate threshold:

| Stage | n_total_epithelial | n_candidate_bipotent | fraction |
|---|---:|---:|---:|
| 9 wpc | 1,799 | 393 | 21.8% |
| 11 wpc | 2,043 | 499 | 24.4% |
| 15 wpc | 1,921 | 335 | 17.4% |
| 18 wpc | 1,812 | 86 | 4.7% |
| 20 wpc | 856 | 2 | 0.2% |
| 22 wpc | 1,086 | 1 | 0.1% |

**Shape:** small early bump 9 → 11 wpc (+2.6 pp), then a strong
strictly monotonic decline 11 → 22 wpc (~250× drop).

The early bump should not be over-interpreted: it is two adjacent
stages, the magnitude is small, and donor-level n is not large
enough to know whether it is biological or sampling.

The post-11-wpc collapse is large and consistent across all three
thresholds and both denominators — that is the load-bearing
observation of v0.

## 5. Threshold robustness

`candidate_bipotent` fraction (% of total epithelial) at all three
SOX2 / SOX9 thresholds — same shape at every threshold:

| Stage | relaxed (>0) | moderate (>0.5) | stringent (>1.0) |
|---|---:|---:|---:|
| 9 wpc | 35.0% | 21.8% | 5.2% |
| 11 wpc | 32.5% | 24.4% | 6.3% |
| 15 wpc | 26.1% | 17.4% | 2.2% |
| 18 wpc | 7.3% | 4.7% | 0.7% |
| 20 wpc | 0.8% | 0.2% | 0.0% |
| 22 wpc | 0.1% | 0.1% | 0.1% |

The **same** 9 → 11 bump and 11 → 22 collapse appears at all three
thresholds. Stricter thresholds compress the absolute fractions
but do not change the trend.

## 6. Denominator robustness

`candidate_bipotent` (moderate threshold) reported two ways. Same
shape under both denominators:

| Stage | as % of total epithelial | as % of candidate_pool |
|---|---:|---:|
| 9 wpc | 21.8% | 28.5% |
| 11 wpc | 24.4% | 35.4% |
| 15 wpc | 17.4% | 23.6% |
| 18 wpc | 4.7% | 6.9% |
| 20 wpc | 0.2% | 2.9% |
| 22 wpc | 0.1% | 1.5% |

The within-pool collapse (28.5% → 1.5%) is on the same trajectory
as the within-total-epithelial collapse, which means the apparent
decline is not just a side-effect of the candidate pool itself
shrinking against committed cells. The candidate pool *does*
shrink, but the bipotent fraction *within* the shrinking pool
collapses too.

(Caveat: at 20 and 22 wpc the within-pool numerator is n=2 and n=1.
The within-pool late-stage percentages are consistent with the
trend but small-n.)

## 7. Curated marker programs

`sc.tl.rank_genes_groups` (Wilcoxon rank-sum) between
`annotated_airway` and `annotated_alveolar` cells, run per stage.
The biological filter in
[`scripts/tipcommit_q1_refinement.py`](../scripts/tipcommit_q1_refinement.py)
flags ribosomal (RPL/RPS), mitochondrial (MT-), and curated
housekeeping rows; 115 of 160 marker rows survive as biological.

DE could only run at 4 of 6 stages — 9 and 11 wpc had **zero**
annotated_alveolar cells.

### 7a. Airway program (recurrent across ≥2 of 4 stages)

Two sub-programs are visible:

- **Mid-stage (15–18 wpc) — secretory + apical/cytoskeletal:**
  AGR2, AGR3 (anterior gradient / secretory; rank 1 every stage),
  CAPS (calcyphosine), HMGN3, S100A11, ANXA1, EZR, DSTN, FXYD3,
  PERP, IFT57. **SOX2 itself recurs at 15 + 18 wpc** (rank 8) — the
  airway TF is recovered as airway-shifted from the canonical
  airway annotation, which is internal validation that the strict
  airway label is biologically coherent.
- **Late-stage (20–22 wpc) — motile-cilia machinery:** CFAP144/276,
  CIMAP1B, CIMIP1, SPMIP6, SPACA9, TSPAN1, DYNLL1, plus
  RSPH1/ZMYND10/DNAH12/CETN2/HOATZ as single-stage hits. Consistent
  with multiciliated-cell maturation.

### 7b. Alveolar program (recurrent across ≥2 of 4 stages)

- **All 4 stages:** **AGER** (canonical AT1), MYL9, CLIC3, SMARCA5
- **3 of 4 stages:** SPARC
- **2 of 4 stages, biologically pointed:** **CLDN18** (alveolar tight
  junction), **SFTPB** (canonical AT2 surfactant), BCAM
- **Single-stage but biologically pointed:** **NKX2-1** at 20 wpc
  (rank 5, alveolar direction); CLDN6, CPM, FOLR1, MDK, AGRN,
  EMP2, S100A10, NREP, SERPINH1 at single stages

Compared with the airway side, the alveolar side has fewer
biological hits per stage — ribosomal genes dominate the alveolar
top-20s at every stage. The biological signal is real (AGER, SFTPB,
CLDN18, NKX2-1) but lower-magnitude than the airway secretory /
ciliated programs.

## 8. Explicit non-claims

The v0 result does **not** claim:

- **Lineage proof.** SOX2 + SOX9 co-expression at scRNA-seq level
  is necessary but not sufficient for true bipotency. Lineage
  tracing / functional validation are wet-lab and out of scope.
  `candidate_bipotent` stays a candidate label throughout.
- **Causation.** The fractions and marker shifts are descriptive.
  No claim that SOX2 or SOX9 *causes* commitment in any given cell.
- **A specific commitment timing.** The fraction declines steeply
  between 11 and 18 wpc, but v0 does not declare "commitment
  becomes irreversible at week N" — that requires lineage tracing.
- **Q2 (gradient vs switch).** The per-stage shape is **monotonic-
  with-an-early-bump**, which is suggestive but not a formal
  resolution of whether the proximal-distal axis is continuous or
  discrete. Q2 is explicitly deferred.
- **Annotation refinement.** v0 trusts the upstream He 2022 cell-
  type labels for committed cells. It does not re-cluster, sub-
  cluster, or correct them.
- **`epithelial cell of lung` = bipotent tip.** This annotation is
  the *candidate* starting pool. The candidate label is preserved
  in every output column.
- **Across-donor statistical inference.** With 10 donors split
  across 6 stages, some stages have 1–2 donors. Per-stage fractions
  are descriptive; no confidence intervals or hypothesis tests.

## 9. Key limitations

1. **No alveolar reference at 9 + 11 wpc.** The dataset has zero
   annotated_alveolar cells at the two earliest stages. Airway-vs-
   alveolar marker DE could not run there. The candidate_bipotent
   fraction itself is still reportable for those stages (no DE
   needed), but the marker comparison only covers 15, 18, 20, 22 wpc.
2. **Tiny candidate pool at 20–22 wpc.** Only 70 (20 wpc) and 66
   (22 wpc) candidate-pool cells. Within-pool fractions collapse
   from numerators of n = 2 and n = 1. The collapse is consistent
   with the trend, but the late-stage within-pool percentages are
   noisy.
3. **20 wpc DE is borderline.** 750 airway vs 36 alveolar cells
   passed the n ≥ 10 threshold for DE, but the imbalance is large.
   20 wpc marker p-values should be read as suggestive, not
   conclusive.
4. **No proximal-distal sub-structure inside the candidate pool.**
   The pool was not sub-clustered. If `epithelial cell of lung`
   contains heterogeneous sub-populations (e.g., proximal-tip vs
   distal-tip), v0 cannot see them.
5. **Single-dataset.** All conclusions rest on one fetal lung atlas.
   No cross-dataset replication (Cao 2020 pan-fetal, LungMAP
   developmental projects, organoid datasets) was attempted.
6. **Donor concentration per stage.** 10 donors across 6 stages
   means some stages are 1–2 donors. Per-stage fractions are
   donor-confounded to an unknown degree.
7. **Two-gene definition.** SOX2 + SOX9 alone defines the bipotent
   class. A multi-gene tip-cell signature (e.g., SOX9 + ID2 + TESC
   + ETV5) would likely be more robust but is v1+ work.

## 10. Where this leaves Q1

**The "when" half of Q1 has a usable v0 answer:** the candidate
SOX2/SOX9 co-expressing fraction in the candidate pool is high
through 11 wpc, declines steeply through 18 wpc, and is essentially
gone by 20–22 wpc — at three independent thresholds and under both
denominator choices. The accompanying marker programs are
biologically coherent (AGR2/3 + cilia for airway; AGER/SFTPB/CLDN18
+ NKX2-1 for alveolar).

**The "how" half (Q2) remains open**, as planned. v0 does not
attempt to resolve whether the transition is gradient-like or
stepwise.
