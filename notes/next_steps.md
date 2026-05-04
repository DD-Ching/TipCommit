# Next Steps — TipCommit

Ordered. Each step has a clear exit condition. Stop after the v0 Q1
MVP per the user's instruction.

## 1. Bootstrap repo scaffold — DONE

- `notes/`, `metadata/`, `scripts/` directories created.
- `notes/evidence_map.md`, `notes/status.md`, `notes/next_steps.md`
  written.
- `git init -b main`; initial commit + push to
  github.com/DD-Ching/TipCommit pending the next bash step.
- Exit condition: scaffold + planning notes visible on GitHub. ✅
  (after the next commit + push).

## 2. Inventory He 2022 in CELLxGENE Census

Goal: determine programmatically whether the He et al. 2022 human
fetal lung dataset is indexed in Census, and what slicing
parameters identify it (collection, dataset_id, donor count, cell
count).

Required outputs:
- A short note `notes/he2022_census_inventory.md` answering:
  - Is He 2022 in the current Census release? Which dataset_id(s)?
  - How many cells, donors, gestational-stage labels?
  - What cell-type ontology terms cover lung epithelial cells
    (CL terms for AT2, AT1, basal, secretory, ciliated, distal tip
    progenitor, etc.)?
  - What development_stage ontology terms appear (HsapDv terms
    for the 5–22 wpc range)?
- A small Python script `scripts/tipcommit_census_inventory.py`
  that runs the lookup and writes the inventory note's data section.

Exit condition: we can answer "what is the Census slice for human
fetal lung epithelial cells in He 2022?" with a concrete dataset_id
and a value_filter expression.

If He 2022 is NOT in the current Census release: fall back to
parallel inventory of Cao 2020 + LungMAP human developmental
projects, and re-decide which is the primary substrate.

## 3. Q1 MVP design — write and lock

Write `notes/q1_mvp_design.md` defining:

- **Census slice:** exact `value_filter` expression and column
  selection.
- **Cell filtering:** include only annotated lung epithelial cells
  (subset by CL terms identified in step 2).
- **Stage binning:** map He 2022's gestational-stage labels into
  3–5 bins (e.g., "early pseudoglandular 5–9 wpc",
  "late pseudoglandular 10–14 wpc", "canalicular 15–18 wpc",
  "saccular 19–22 wpc"; final bins TBD by what data exists).
- **Expression scoring:** log-normalised SOX2 and SOX9 expression
  per cell. Define SOX2+ / SOX9+ thresholds based on the
  observed distribution within epithelial cells (e.g., > median
  + 1 MAD, or > 0 in log space — to be locked after a histogram
  inspection).
- **Cell classification:** per cell, assign one of:
  - `bipotent_double_positive` (SOX2+ SOX9+)
  - `airway_committed` (SOX2+ SOX9−)
  - `alveolar_committed` (SOX2− SOX9+)
  - `double_negative` (other / non-progenitor or AT1-like)
- **Per-stage output:** fraction of cells in each class at each
  stage.
- **Marker comparison:** for each stage, run
  `sc.tl.rank_genes_groups` between `airway_committed` and
  `alveolar_committed` cells; report the top 10–20 genes per
  direction.
- **Output files:**
  - `metadata/q1_per_stage_fractions.csv` — wide format
    (stage × class).
  - `metadata/q1_long_format.csv` — figure-ready (stage, class, count, fraction).
  - `metadata/q1_markers_by_stage.csv` — stage × direction × top-N
    genes with log-fold-change and adjusted p-value.

Exit condition: design note committed; reviewer (the user) can
confirm the MVP shape before implementation begins.

## 4. Implement Q1 MVP

Single new script `scripts/tipcommit_q1_mvp.py`:

- Reads from Census (one slice per run; cache locally if the slice
  is large).
- Applies the design's filtering + scoring + classification.
- Writes the three CSVs.
- Prints a console summary (stage × class fractions table).

Dependencies installed at this step: `cellxgene-census`, `scanpy`,
`anndata`, plus their transitive deps (numpy, pandas, scipy). No
other packages.

Exit condition: script runs end-to-end; three CSVs in `metadata/`;
console summary shows expected magnitudes (e.g., bipotent fraction
declines with stage).

## 5. Stop and report

After step 4, the user has the Q1 MVP. Stop. Report:

- The headline numbers: bipotent fraction at earliest vs latest
  stage.
- Whether the per-stage trend is monotonic (suggesting gradual
  commitment) or step-wise (suggesting a discrete switch — but Q2
  isn't being formally tested in v0).
- Top marker-gene shifts by stage.
- What v0 explicitly does *not* claim (commitment timing in any
  particular cell, statistical confidence, functional validation).

## Constraints carried forward

- **Q1 only.** Q2 is deferred. Do not run pseudotime, trajectory
  inference, or any axis-continuity test in v0.
- **No broad platform.** This is one Q1-focused script, not a
  general single-cell-analysis pipeline.
- **No packages beyond the cellxgene-census + scanpy + anndata
  stack.** No scvi-tools, no plotnine, no seaborn. Use
  matplotlib only if a figure is essential to v0 (probably not —
  the figure-ready table is enough).
- **Plan before code.** Each new analysis step gets a 1–2 paragraph
  note in `notes/` before implementation.
- **Small commits, push after each meaningful step.** Same cadence
  as Gain.
- **Separate known biology / inference / hypothesis.** Every note
  must label which is which.
