# Q1.5b Replication Design

**Date:** 2026-05-05
**Status:** design contract for the replication pass. No
implementation in this turn — design first, implement after.
**Predecessor:** Q1 v0 ([`notes/q1_summary.md`](q1_summary.md)).

This is a **single replication pass on one independent dataset**.
The biological scope is unchanged. The method is unchanged. Only the
substrate changes.

---

## A. Chosen replication dataset

**Cao et al. 2020 — "A human cell atlas of fetal gene expression"
(*Science*, doi:10.1126/science.aba7721)** — accessed via CELLxGENE
Census, dataset_id **`fa27492b-82ff-4ab7-ac61-0e2b184eee67`**
(the curated "1 million cells subset" Census release).

| Attribute | Value |
|---|---|
| Census dataset_id | `fa27492b-82ff-4ab7-ac61-0e2b184eee67` |
| Title | "Survey of human embryonic development (1 million cells subset)" |
| Collection | "A human cell atlas of fetal gene expression" |
| Publication | doi:10.1126/science.aba7721 (Cao et al. 2020) |
| Lung-subset cells (`tissue_general == 'lung'`) | 53,429 |
| Donors | 11 |
| Assay | sci-RNA-seq3 (combinatorial barcoding, single-nucleus-style) |
| Developmental stages present in lung | 12, 13, 15, 16, 17 wpc |

### Why this dataset is the best replication substrate

1. **Independent of He 2022.** Different lab (Shendure/Cao vs
   Quake/He), different assay platform (sci-RNA-seq3 vs 10x 5′ v1),
   different annotation pipeline. A trend that survives both
   substrates is not an artefact of either single pipeline.
2. **Pan-organ atlas, not lung-targeted.** Cao 2020 was annotated
   in a multi-organ context — its lung cells weren't pre-selected
   for proximal-distal-axis biology. That removes a possible
   selection bias from the He 2022 collection design.
3. **Same Census interface.** The `cellxgene-census` value_filter
   pattern from the Q1 MVP works directly: only the dataset_id
   string changes. No new dependencies, no new data formats.
4. **Tractable size.** ~53k lung cells fetches in 1–2 minutes; the
   full 4M-cell parent atlas (`f7c1c579-2dc0-47e2-ba19-8165c5a0e353`,
   ~218k lung cells) is the canonical Cao 2020 release but is
   ~4× heavier and the 1M subset is the same scientific data
   downsampled. Per the user's "smallest replication pass" rule,
   the 1M subset wins.

### Why other candidates were rejected

- **He 2022 organoid sister-dataset (`4023a2bc...`).** Same
  collection — not independent of the primary substrate. Useful as
  a v1+ in vitro cross-check but not a replication.
- **He 2022 sub-cell-type splits (Epithelium-no-cilium / Cilium /
  Myeloid / Endothelium / etc.).** All sub-slices of the He 2022
  "All cells" dataset. Not independent.
- **"Early human lung immune cell development" datasets
  (`350237e0`, `fcadb222`).** Both immune-focused (the larger has
  670k cells but only ~869 epithelial-keyword hits per the Step 2
  inventory). Wrong cell-population focus.
- **LungMAP "Single-cell multiomic profiling of human lungs across
  age groups" (`3de0ad6d`).** Only 31 wpc fetal stage represented
  — out of the 9–22 wpc Q1 window. Wrong stage range.
- **Han 2020 "Human Cell Landscape" (`2adb1f8a`).** Only 11 + 12
  wpc lung stages and only ~9.6k lung cells. Stage range too narrow
  for a trend test.

### Why NOT use the full Cao 2020 atlas (`f7c1c579`)

- 4M-cell parent has ~218k lung cells (vs 53k in the subset) — 4×
  longer fetch, no biological gain for the replication question.
- The 1M subset is `is_primary_data=False` in Census because it
  duplicates cells in the parent — a Census-internal flag for
  deduplication in aggregations, not a data-quality flag. The cells
  are real primary measurements from the same Cao 2020 paper.
- If the replication is inconclusive at 53k cells, the full atlas
  is the v1 follow-up; the user can choose to re-run on `f7c1c579`
  by changing one constant.

---

## B. How closely the original Q1 definitions can be reused

Cao 2020's lung-subset epithelial annotation is **much coarser**
than He 2022's. Only 4 epithelial cell types appear:

| Cao 2020 cell type | n cells (1M subset, lung) | Closest He 2022 analog |
|---|---:|---|
| `epithelial cell of lower respiratory tract` | 24,809 | `epithelial cell of lung` (the candidate pool) |
| `ciliated epithelial cell` | 659 | `lung multiciliated epithelial cell` (annotated airway) |
| `squamous epithelial cell` | 70 | `squamous epithelial cell` (annotated airway, outlier) |
| `neuroendocrine cell` | 346 | (no direct analog; excluded) |

**No** AT1, AT2, basal, club, or lung secretory annotations exist in
Cao 2020's lung subset. The annotation pipeline simply did not
sub-divide the alveolar / airway lineages at this resolution.

### Reuse map — Q1 v0 (He 2022) → Q1.5b (Cao 2020)

| Q1 v0 class | Q1.5b analog | Reuse status |
|---|---|---|
| `candidate_pool` = `epithelial cell of lung` | `epithelial cell of lower respiratory tract` | **direct analog** |
| `candidate_bipotent` (SOX2 > T ∧ SOX9 > T within pool) | same rule | **direct** |
| `transitioning_airway` (SOX2 > T only within pool) | same rule | **direct** |
| `transitioning_alveolar` (SOX9 > T only within pool) | same rule | **direct** |
| `undefined` (neither within pool) | same rule | **direct** |
| `annotated_airway` = basal + club + secretory + multiciliated + respiratory + squamous | `ciliated epithelial cell` + `squamous epithelial cell` only | **partial** — only ciliated+squamous; no basal/club/secretory annotated |
| `annotated_alveolar` = AT1 + AT2 | (none) | **CANNOT REPLICATE** |

### Threshold sweep — direct reuse

Same three thresholds (relaxed > 0, **moderate > 0.5 default**,
stringent > 1.0). Same normalisation
(`sc.pp.normalize_total(target_sum=1e4)` + `sc.pp.log1p`). Same
candidate-not-proven discipline applied to `candidate_bipotent`.

### Marker DE — must change

The He 2022 MVP ran Wilcoxon between `annotated_airway` and
`annotated_alveolar`. With **zero** annotated alveolar cells in
Cao 2020, that exact comparison is impossible. Two options were
considered:

- **Option A**: skip marker DE entirely; report fractions only.
- **Option B (chosen)**: run Wilcoxon between
  `transitioning_airway` (SOX2 > 0.5 only, within candidate pool)
  and `transitioning_alveolar` (SOX9 > 0.5 only, within candidate
  pool), at each stage. Both groups are defined within the same
  candidate pool by the same expression rule.

**Caveat carried forward**: SOX2 and SOX9 themselves trivially
separate these two groups by construction. They will be excluded
from the curated marker list in the results note. Other genes
(AGR2/3, CAPS, AGER, SFTPB, CLDN18, NKX2-1) are not part of the
defining rule, so their recovery is informative.

This Option B comparison is **not** the same as the He 2022 strict
annotation-based comparison. The replication on the marker side is
*directional* (do airway-leaning vs alveolar-leaning cells split by
similar gene programs?) rather than *exact* (do the same top-N gene
lists appear?).

---

## C. What must change because of annotation differences (summary)

1. **No `annotated_alveolar` group.** Drop from outputs; flag as
   "not annotated in Cao 2020" wherever the schema would expect it.
2. **`annotated_airway` shrinks** to ciliated + squamous only
   (~729 cells in the lung subset vs ~3,257 in He 2022's airway
   set). Per-stage counts will be sparse.
3. **Marker DE switches to transitioning_airway vs
   transitioning_alveolar** within the candidate pool, with the
   SOX2/SOX9 circularity caveat above.
4. **Stage coverage shrinks** from He's 9–22 wpc (6 stages) to
   Cao's 12–17 wpc (5 stages). Direct stage overlap with He 2022:
   only **15 wpc**.
5. **Assay difference**: sci-RNA-seq3 (Cao) vs 10x 5′ v1 (He).
   Both are short-read scRNA-seq with similar dynamic range, but
   sci-RNA-seq3 typically has lower per-cell UMI counts → higher
   dropout → potentially fewer SOX2/SOX9 detection events at the
   same threshold. The threshold sweep handles this; the moderate
   default (>0.5) may catch fewer cells in Cao than in He at the
   same stage. This is expected and not a failure mode.

---

## D. Successful-replication criteria

A **successful replication** requires all three of:

1. **Direction.** Within the Cao 2020 candidate pool, the
   `candidate_bipotent` fraction at the moderate threshold shows
   a **declining** trend across stages (12 → 17 wpc). Strict
   monotonicity is not required; the trend should be visibly
   downward, not flat or rising.
2. **Magnitude (rough).** At 15 wpc — the only direct stage
   overlap — Cao 2020's within-pool `candidate_bipotent` fraction
   sits in the **same order of magnitude** as He 2022's 23.6%
   (acceptable range: roughly 5%–60% — a much lower assay-driven
   value is fine; a much higher one would suggest something
   different).
3. **Marker direction.** At least 2 of the curated airway markers
   (AGR2, AGR3, CAPS) appear in the top-20 transitioning_airway
   list **or** at least 2 of the curated alveolar markers (AGER,
   SFTPB, CLDN18) appear in the top-20 transitioning_alveolar
   list, at one or more stages.

If all three hold → **replicates.**

## E. Non-replication criteria

A **failed replication** is any one of:

- `candidate_bipotent` fraction is **flat** (≤ 2-fold range across
  stages) at the moderate threshold within the candidate pool.
- `candidate_bipotent` fraction **rises** with stage (opposite
  trend).
- 15 wpc Cao value is order-of-magnitude different from He 2022
  (e.g., < 2% or > 70%) AND no other Cao stage is in range either.
- None of {AGR2, AGR3, CAPS, AGER, SFTPB, CLDN18} appear in the
  Cao top-20 lists at any stage.

A **partial replication** is direction-replicated (criterion 1)
but missing magnitude (criterion 2) or marker direction (criterion
3). Partial replication still strengthens the v0 finding but with
appropriate caveats.

---

## F. What this replication explicitly does NOT test

- **The 9–11 wpc early bump.** Cao 2020 lung doesn't cover those
  stages. The bump remains a single-dataset observation.
- **The 18–22 wpc late collapse.** Cao 2020 stops at 17 wpc.
  Whether the bipotent fraction goes to ~0.1% by 22 wpc remains a
  single-dataset observation.
- **The strict annotated_airway vs annotated_alveolar marker
  programs.** Cao 2020 doesn't have annotated alveolar cells; the
  marker comparison is on a different (expression-defined) split.
- **Whether `epithelial cell of lower respiratory tract` =
  `epithelial cell of lung`.** The two annotations come from
  different ontology / annotation pipelines. The replication
  treats them as analogs but doesn't prove cell-by-cell equivalence.
- **Donor / batch effects.** No donor-stratified analysis. With
  11 Cao donors across 5 stages, donor confounding mirrors the
  He 2022 limitation.

---

## G. Reuse of TipCommit logic

The replication script will be a near-clone of
`scripts/tipcommit_q1_mvp.py` with the constants below changed.
**No new functions, no new method families, no broadening.**

Constant changes:

```python
DATASET_ID            = "fa27492b-82ff-4ab7-ac61-0e2b184eee67"  # Cao 2020 1M
CANDIDATE_POOL        = "epithelial cell of lower respiratory tract"
ANNOTATED_AIRWAY      = {"ciliated epithelial cell", "squamous epithelial cell"}
ANNOTATED_ALVEOLAR    = set()  # absent in Cao 2020 lung
EPITHELIAL_CELL_TYPES = ANNOTATED_AIRWAY | {CANDIDATE_POOL}  # neuroendocrine excluded
MARKER_GROUP_A        = "transitioning_airway"        # was annotated_airway
MARKER_GROUP_B        = "transitioning_alveolar"      # was annotated_alveolar
```

Logic changes:

- The classify-cells function still emits the same six labels;
  `annotated_alveolar` will simply be empty.
- The per-stage CSV will keep all the same columns; alveolar-side
  counts will be 0 throughout.
- The markers function will compare the two transitioning-* groups
  instead of the two annotated-* groups, with the methodological
  caveat noted above.
- Curation of markers (gene_class flagging) reuses the
  `tipcommit_q1_refinement.py` logic by shared constants /
  inline copy. No new file.

Output filenames are namespaced with `q15b_replication_*` so they
sit alongside the v0 outputs without overwriting them.

---

## H. Implementation outline (for the next turn)

A single new script `scripts/tipcommit_q15b_replication.py` that:

1. Opens Census, pulls AnnData per the new value_filter.
2. Filters to the 4 epithelial Cao annotations (excluding
   neuroendocrine).
3. Normalises (`sc.pp.normalize_total(target_sum=1e4)` + `sc.pp.log1p`).
4. Computes per-cell SOX2 / SOX9 log-normalised expression, classes
   at the three thresholds.
5. Per-stage aggregation; writes
   `metadata/q15b_replication_per_stage_fractions.csv` (wide) and
   `metadata/q15b_replication_long_format.csv` (long, with
   threshold label).
6. Per stage, runs Wilcoxon between `transitioning_airway` and
   `transitioning_alveolar` cells (defined at the moderate threshold);
   writes `metadata/q15b_replication_markers_by_stage.csv` with
   `gene_class` already attached (no separate refinement script).
7. Prints fraction tables + curated marker top hits.

Estimated script length: ~280 LoC. Estimated runtime: 1–3 min
(Census fetch dominates).

After the script writes the three CSVs, write
`notes/q15b_replication.md` with the actual numbers and a verdict
against the criteria in section D / E.

No new packages, no plotting, no Q2.
