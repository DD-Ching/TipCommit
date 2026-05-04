# Q1 MVP Design

**Status:** design contract for Step 4 (Q1 MVP implementation). No
code in this turn — Step 3 only per `notes/next_steps.md`.

**Question:** when and how does the bipotent SOX2/SOX9 co-expressing
tip cell commit to airway vs alveolar fate in human fetal lung?

**Substrate (locked in Step 2):** He 2022 atlas in CELLxGENE Census,
dataset `3dc61ca1-ce40-46b6-8337-f27260fd9a03`. See
[`notes/he2022_census_inventory.md`](he2022_census_inventory.md).

**User correction carried into this design:** "epithelial cell of
lung" is **a candidate starting pool**, NOT proven bipotent. The
design treats this annotation as the substrate within which SOX2 +
SOX9 expression-based sub-classification identifies *candidate*
bipotent cells. The MVP does not claim bipotency; it reports
expression-based candidates and their downstream marker behaviour.

---

## A. Substrate

### Exact Census slice

```python
value_filter = (
    "tissue_general == 'lung' and "
    "dataset_id == '3dc61ca1-ce40-46b6-8337-f27260fd9a03'"
)
```

Pulled with `cellxgene_census.get_anndata(census, organism="Homo sapiens", obs_value_filter=value_filter)`.
Returns an AnnData object with all 62,759 fetal-stage cells (single
assay 10x 5′ v1, 100% primary data).

### Exact stage labels — all 6, unbinned

Keep all six developmental_stage labels at their native granularity
(week-by-week). The cell counts per stage are large enough to
support per-stage fraction estimates without binning, and binning
would lose the timing information that Q1 most cares about.

| Stage | Total cells | Approx. epithelial (15%) |
|---|---:|---:|
| 9th wpc | 8,280 | ~ 1,240 |
| 11th wpc | 9,531 | ~ 1,430 |
| 15th wpc | 12,071 | ~ 1,810 |
| 18th wpc | 13,170 | ~ 1,975 |
| 20th wpc | 7,490 | ~ 1,125 |
| 22nd wpc | 12,217 | ~ 1,830 |

### Exact epithelial populations to include / exclude

**Include** (9 cell-type annotations from He 2022; total ~9,517 cells):

| cell_type | Cells | Treatment in Q1 |
|---|---:|---|
| `epithelial cell of lung` | 5,595 | **Candidate starting pool** — sub-classified by SOX2/SOX9 expression |
| `lung multiciliated epithelial cell` | 1,641 | annotated airway (committed) |
| `respiratory tract epithelial cell` | 991 | annotated airway (committed; broad annotation) |
| `pulmonary alveolar type 2 cell` | 382 | annotated alveolar (committed; AT2) |
| `pulmonary alveolar type 1 cell` | 283 | annotated alveolar (committed; AT1) |
| `basal cell` | 278 | annotated airway (committed; basal) |
| `club cell` | 208 | annotated airway (committed; club) |
| `lung secretory cell` | 128 | annotated airway (committed; secretory) |
| `squamous epithelial cell` | 11 | annotated airway (committed; outlier — included for completeness) |

**Exclude:**
- All non-epithelial cell types (fibroblasts, smooth muscle, endothelial,
  immune, neuroendocrine, etc.). Q1 is strictly about epithelial-pool
  commitment.

---

## B. Working cell-state definitions

Hybrid annotation + expression rule. Six classes per cell.

### Annotation-based (canonical committed labels — trust the upstream paper)

- **`annotated_airway`** = cell whose `cell_type` is one of:
  `lung multiciliated epithelial cell`, `respiratory tract epithelial cell`,
  `basal cell`, `club cell`, `lung secretory cell`,
  `squamous epithelial cell`. Total ≈ 3,257 across the dataset.
- **`annotated_alveolar`** = cell whose `cell_type` is one of:
  `pulmonary alveolar type 1 cell`, `pulmonary alveolar type 2 cell`.
  Total ≈ 665.

These are the most committed cells in the dataset. Their annotations
came from the upstream paper's clustering and marker analysis; the
MVP trusts those labels.

### Expression-based (sub-classification of the candidate pool)

Within the **`epithelial cell of lung`** annotation (≈ 5,595 cells —
the candidate starting pool):

- **`candidate_bipotent`** = SOX2+ AND SOX9+. The SOX2/SOX9 co-
  expressing subset of the candidate pool. **Candidate** because
  co-expression at scRNA-seq level is necessary but not sufficient
  for true bipotency; functional / lineage-tracing validation is
  out of scope.
- **`transitioning_airway`** = SOX2+ AND SOX9−. Within the candidate
  pool, cells leaning proximal/airway by expression but not yet
  annotated as a committed airway cell-type by the upstream paper.
- **`transitioning_alveolar`** = SOX2− AND SOX9+. Within the
  candidate pool, cells leaning distal/alveolar by expression but
  not yet annotated as AT1 or AT2.
- **`undefined`** = SOX2− AND SOX9−. Within the candidate pool,
  cells that don't express either marker above threshold. Could be
  later-progenitor states, dropouts, or non-progenitor noise.

### The user's three named groups, mapped to these classes

The user's "candidate SOX2+/SOX9+ co-expressing pool" =
**`candidate_bipotent`** above.

The user's "airway-leaning committed cells" can be reported two
ways in the MVP outputs:

- **strict_airway_committed** = `annotated_airway` only (canonical
  committed cells from the upstream annotation).
- **broad_airway_leaning** = `annotated_airway` ∪ `transitioning_airway`
  (canonical + expression-leaning).

The user's "alveolar-leaning committed cells":

- **strict_alveolar_committed** = `annotated_alveolar` only.
- **broad_alveolar_leaning** = `annotated_alveolar` ∪ `transitioning_alveolar`.

The MVP CSVs report both strict and broad variants so the user can
inspect either reading.

---

## C. Threshold choices

### Primary thresholds (locked for v0)

After standard Scanpy normalisation (`sc.pp.normalize_total(target_sum=1e4)`
followed by `sc.pp.log1p`), per-cell expression is in log-normalised
units.

- **SOX2+** ≡ log-normalised `SOX2` expression **> 0** (i.e., any
  non-zero count survived normalisation).
- **SOX9+** ≡ log-normalised `SOX9` expression **> 0**.
- **Co-expression (`candidate_bipotent`)** ≡ both SOX2+ AND SOX9+
  per the rules above.

This is the **most permissive** threshold — any cell with at least
one detected transcript of each gene crosses it. Justified for v0
because:

1. scRNA-seq dropouts are heavy; stricter thresholds will under-
   count co-expressing cells.
2. The rule is fully explainable and reproducible.
3. The upstream paper's annotation work has already validated
   these cells as epithelial; we're not trying to detect epithelial
   identity from SOX2/SOX9 alone.

### Sensitivity checks (run alongside the primary)

Re-run the per-stage fractions at two stricter thresholds:

- **`threshold_relaxed`** ≡ > 0 (the primary; documented above)
- **`threshold_moderate`** ≡ > 0.5 (≈ 1.6 normalised counts)
- **`threshold_stringent`** ≡ > 1.0 (≈ 2.7 normalised counts)

The MVP reports `candidate_bipotent` fraction at all three
thresholds in the per-stage CSV. If the **trend across stages** is
robust to threshold choice (e.g., bipotent fraction declines with
stage at all three thresholds), that is a robust observation. If the
trend disappears at stricter thresholds, the v0 finding is
threshold-driven and must be flagged.

### Distribution-aware sanity check (logged, not parameterised)

Before locking the threshold for the per-stage analysis, the script
will print SOX2 and SOX9 expression-distribution summaries within
the candidate pool (5-number summary + fraction at each threshold).
This is a sanity check, not a parameter — the design's locked
thresholds (> 0, > 0.5, > 1.0) hold regardless of what the
distribution looks like, but seeing the distribution lets the user
spot if e.g. SOX2 is bimodal (which would suggest a different
threshold is more natural for v1).

---

## D. MVP outputs

Three CSVs. Schemas locked here.

### `metadata/q1_per_stage_fractions.csv` (wide format)

One row per stage. Columns:

| column | type | description |
|---|---|---|
| `stage` | string | e.g. `"15th week post-fertilization stage"` |
| `n_total_epithelial` | int | all 9 epithelial cell types at this stage |
| `n_annotated_airway` | int | strict airway committed |
| `n_annotated_alveolar` | int | strict alveolar committed |
| `n_candidate_pool_total` | int | "epithelial cell of lung" annotation count |
| `n_candidate_bipotent_relaxed` | int | SOX2+ ∧ SOX9+ at threshold > 0 |
| `n_candidate_bipotent_moderate` | int | SOX2+ ∧ SOX9+ at threshold > 0.5 |
| `n_candidate_bipotent_stringent` | int | SOX2+ ∧ SOX9+ at threshold > 1.0 |
| `n_transitioning_airway` | int | candidate pool ∩ SOX2+SOX9− at relaxed threshold |
| `n_transitioning_alveolar` | int | candidate pool ∩ SOX2−SOX9+ at relaxed threshold |
| `n_undefined` | int | candidate pool ∩ SOX2−SOX9− at relaxed threshold |
| `frac_<each_class>` | float | n / n_total_epithelial |

6 rows. One file. Wide format, easy to inspect.

### `metadata/q1_long_format.csv` (figure-ready)

Long format; one row per (stage, class) combination. Columns:

| column | type | description |
|---|---|---|
| `stage` | string | as above |
| `class` | string | one of: `annotated_airway`, `annotated_alveolar`, `candidate_bipotent_relaxed`, `transitioning_airway`, `transitioning_alveolar`, `undefined` |
| `n_cells` | int | count |
| `fraction` | float | n_cells / n_total_epithelial at this stage |

6 stages × 6 classes = 36 rows. Used by any downstream plotting
tool (matplotlib, seaborn, R ggplot2, etc.) — the MVP does not
ship plots itself.

### `metadata/q1_markers_by_stage.csv` (marker shifts)

Top genes differentiating `annotated_airway` vs `annotated_alveolar`
at each stage, via `sc.tl.rank_genes_groups` (Wilcoxon rank-sum,
which is Scanpy's default for this method). Columns:

| column | type | description |
|---|---|---|
| `stage` | string | as above |
| `direction` | string | `airway_top` (genes up in airway vs alveolar) or `alveolar_top` (up in alveolar vs airway) |
| `rank` | int | 1 to 20 (top 20 per direction per stage) |
| `gene` | string | gene symbol |
| `log_fold_change` | float | log2 fold-change reported by Scanpy |
| `adjusted_p_value` | float | Benjamini-Hochberg-adjusted p-value |

6 stages × 2 directions × top 20 = up to 240 rows. Stages with
fewer than ~20 cells in either group are skipped (with a logged
warning).

---

## E. Minimal figure-ready result

The MVP must produce **one core observation** plus **two supporting
ones** that a reader can extract from the CSVs without re-running
the script.

### Core observation

**Per-stage fraction of candidate SOX2+/SOX9+ cells** (the
`candidate_bipotent_relaxed` class as a fraction of `n_total_epithelial`),
plotted across the 6 stages (9 → 22 wpc).

The expected shape (a hypothesis to be reported as such, not a
finding): a monotonic decline as the bipotent pool resolves into
committed populations. If the actual shape differs (e.g., spike in
the middle, no decline), that is itself a finding.

### Supporting observations

1. **Per-stage airway-leaning vs alveolar-leaning fractions** — both
   strict (`annotated_airway`, `annotated_alveolar`) and broad
   (annotated ∪ transitioning) variants. Should show a complementary
   rise as the bipotent fraction declines.
2. **Top marker shifts associated with each direction** — at each
   stage, the top 20 airway-up vs top 20 alveolar-up genes from
   `sc.tl.rank_genes_groups`. Stage-to-stage stability of these
   markers is itself informative (stable markers = the labels are
   coherent across stages; shifting markers = stage-specific biology
   or annotation drift).

### What the MVP does NOT plot

The MVP ships the figure-ready long-format CSV. It does **not**
build matplotlib / seaborn / plotnine plotting infrastructure. The
user can plot from `q1_long_format.csv` in any tool. This is per
the user's "do not build plotting infrastructure" rule.

---

## F. Explicit non-claims

The MVP will NOT claim:

- **No proof of true lineage commitment.** Co-expression of SOX2 and
  SOX9 at scRNA-seq level is necessary but not sufficient for
  bipotency. Functional / lineage-tracing validation is wet-lab
  territory and out of scope.
- **No causal inference.** The fractions and marker shifts are
  descriptive. No claim that SOX2 or SOX9 *causes* commitment in any
  particular cell.
- **No Q2 gradient-vs-switch claim.** Whether the proximal-distal
  axis is continuous or discrete is Q2; explicitly deferred. Even if
  the per-stage fractions look monotonic-and-gradual or step-wise-
  and-abrupt, the MVP will not classify the axis as continuous or
  discrete.
- **No pseudotime-heavy interpretation.** No `sc.tl.dpt`,
  `scanpy.external.tl.palantir`, scvelo, or similar. The MVP's only
  ordering is the dataset's `development_stage` label (donor age in
  wpc). If trajectory inference is justified later, it gets its own
  design note.
- **No annotation refinement.** The MVP trusts the He 2022 paper's
  cell-type labels for committed cells (basal, club, AT1, AT2,
  etc.). It does not re-cluster, re-annotate, or correct the
  upstream annotation.
- **No claim that "epithelial cell of lung" = bipotent tip.**
  This annotation is the **candidate starting pool**; SOX2/SOX9
  expression sub-classifies it into candidate bipotent / transitioning /
  undefined. The candidate label is preserved throughout the
  outputs.
- **No claim about specific commitment timing.** Even if the
  candidate bipotent fraction declines monotonically, the MVP
  reports the trend; it does not declare "commitment becomes
  irreversible at week N" — that requires lineage tracing.
- **No across-donor statistical inference.** With 10 donors split
  across 6 stages, some stages have 1–2 donors. Per-stage fractions
  are descriptive; no confidence intervals or hypothesis tests in
  v0.

---

## Implementation outline (for Step 4 — NOT this turn)

A single script `scripts/tipcommit_q1_mvp.py` that:

1. Opens Census, pulls AnnData per the locked `value_filter`.
2. Filters to the 9 epithelial cell types.
3. Normalises (`sc.pp.normalize_total(target_sum=1e4)` + `sc.pp.log1p`).
4. Computes per-cell SOX2 and SOX9 log-normalised expression and
   thresholds at `> 0`, `> 0.5`, `> 1.0`.
5. Assigns each cell to one of 6 classes per the rules in section B.
6. Aggregates per stage; writes `q1_per_stage_fractions.csv` (wide)
   and `q1_long_format.csv` (long).
7. Per stage, runs `sc.tl.rank_genes_groups` between
   `annotated_airway` and `annotated_alveolar` cells; writes
   `q1_markers_by_stage.csv`.
8. Prints the core observation table to stdout.

Estimated script length: 200–300 LoC stdlib + the cellxgene-census /
scanpy / anndata stack. Estimated runtime: 1–3 min (Census fetch is
the main cost).

No new packages beyond the locked stack. No plotting code.

## Out of scope for v0

- Q2 (axis continuity) — separate design when Q1 ships.
- Pseudotime / trajectory inference / RNA velocity.
- Cell-type re-annotation or sub-clustering of "epithelial cell of
  lung".
- Cross-dataset integration (other fetal lung atlases).
- Mouse comparison.
- Functional validation links.
- Plotting infrastructure.
- scvi-tools.
