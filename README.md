# TipCommit

> **Q1 result on bipotent SOX2 / SOX9 commitment in human fetal lung
> epithelium.** Census + Scanpy. Q1 v0 cycle complete.

A small expression-correlation project asking the **"when"** half of
a single biological question:

> When and how does the bipotent SOX2 / SOX9 co-expressing tip cell
> commit to airway vs alveolar fate in **human** fetal lung?

TipCommit is the successor to
[Gain v0.1.0](https://github.com/DD-Ching/Gain/releases/tag/v0.1.0)
(an evidence-audit / substrate-gap repo whose chain stopped at the
public-data ChIP gap). TipCommit picks up the
expression-correlation route documented in
[Gain's next-project decision memo](https://github.com/DD-Ching/Gain/blob/main/notes/next_project_decision.md),
scoped strictly to Q1.

## Project status

**Q1 v0 cycle complete (2026-05-05). Standing outputs frozen.**
No further Q1 refinement passes are planned in this repo phase.
Q2 (continuous gradient vs discrete switch) is explicitly
deferred — see *What remains open* below.

### Start here (in order)

1. **[`notes/q1_key_findings.md`](notes/q1_key_findings.md)** —
   eight bullet points, two-minute read.
2. **[`notes/q1_summary.md`](notes/q1_summary.md)** — the full v0
   narrative for a new researcher (question, dataset, definitions,
   per-stage trend, threshold + denominator robustness, curated
   marker programs, non-claims, limitations).
3. **[`notes/he2022_census_inventory.md`](notes/he2022_census_inventory.md)** —
   substrate inventory locking the He 2022 atlas as the Q1 dataset.
4. **[`notes/q1_mvp_design.md`](notes/q1_mvp_design.md)** — design
   contract for the MVP (cell-state definitions, threshold sweep,
   output schemas, non-claims).
5. **[`notes/q1_refinement.md`](notes/q1_refinement.md)** — Q1.1
   refinement pass (cleaner state breakdown, curated markers,
   trend-robustness check).

## What this repo establishes

1. **A robust per-stage trend in human fetal lung.** The candidate
   SOX2 / SOX9 co-expressing fraction inside the He 2022 atlas's
   candidate epithelial pool collapses from ~24% at 11 wpc to ~0.1%
   by 22 wpc at the moderate threshold. The trajectory is the same
   under three independent thresholds (relaxed / moderate /
   stringent) and under both denominator choices (total epithelial
   vs candidate pool only). This is the load-bearing finding.
2. **Coherent stage-structured marker programs.** Airway-shifted
   markers organise into a 15–18 wpc secretory / apical program
   (AGR2, AGR3, CAPS, ANXA1, EZR, PERP) and a 20–22 wpc
   motile-cilia program (CFAP144/276, CIMAP1B, CIMIP1, TSPAN1,
   DYNLL1), with SOX2 itself recurring at 15 + 18 wpc. Alveolar-
   shifted markers recover the canonical AT1 / AT2 set (AGER, SFTPB,
   CLDN18) plus NKX2-1 at 20 wpc.
3. **Honest small-n boundaries.** Two earliest stages (9 + 11 wpc)
   have **zero** annotated alveolar cells in this dataset, so the
   marker DE only runs at 4 of 6 stages. Two latest stages (20 + 22
   wpc) have very small candidate pools (n = 70, 66). Both
   limitations are documented in every output table that touches
   those stages.
4. **A reproducible MVP + refinement pipeline.** Two scripts
   ([`scripts/tipcommit_q1_mvp.py`](scripts/tipcommit_q1_mvp.py),
   [`scripts/tipcommit_q1_refinement.py`](scripts/tipcommit_q1_refinement.py))
   that run end-to-end against a public Census release on the
   `cellxgene-census` + `scanpy` + `anndata` stack with no other
   dependencies.

## What this repo is *not*

- **Not a lineage-tracing claim.** SOX2 + SOX9 co-expression at
  scRNA-seq level is necessary but not sufficient for true
  bipotency. `candidate_bipotent` stays a candidate label
  throughout.
- **Not a Q2 result.** v0 does not test whether the proximal-distal
  axis is continuous or discrete.
- **Not a multi-dataset replication.** Conclusions rest on one
  fetal lung atlas (He 2022). No cross-dataset re-derivation.
- **Not a perturbation, spatial, or pseudotime analysis.** The
  ordering is the dataset's `development_stage` label.

## What remains open

- **The "how" half of Q1 (= Q2).** Whether the bipotent → committed
  transition is a continuous gradient or a discrete switch. The v0
  per-stage shape is *monotonic-with-an-early-bump*, which is
  suggestive but not a formal answer.
- **Sub-structure inside the candidate pool.** The
  `epithelial cell of lung` annotation was treated as one pool; if
  it contains proximal-vs-distal-tip heterogeneity, v0 cannot see
  it.
- **Cross-dataset replication.** Cao 2020 pan-fetal, LungMAP
  developmental projects, and the He 2022 organoid sister-dataset
  could each independently test the per-stage decline.
- **Multi-gene tip signatures.** A SOX9 + ID2 + TESC + ETV5
  composite would likely be more robust than the two-gene SOX2 +
  SOX9 rule used here.
- **Donor-level resampling / per-donor variance.** With 10 donors
  across 6 stages, some stages are 1–2 donors. v0 reports
  fractions descriptively only.

## Recommended next question

The most natural next step is **Q2** (axis continuity), but two
softer Q1.5-style options sit between Q1 and Q2:

- **Q1.5a — sub-structure of the candidate pool.** Re-cluster
  `epithelial cell of lung` cells (5,943 in this run) and ask
  whether the candidate_bipotent subset corresponds to a single
  coherent sub-cluster or distributes across several.
- **Q1.5b — cross-dataset replication.** Re-run the same pipeline
  on a second fetal lung dataset (organoid or Cao 2020 lung subset)
  and check whether the per-stage decline reproduces.

Either of these would tighten the v0 result before committing to
the full Q2 design. The decision is left for the next session.

## Substrate

| Resource | Use |
|---|---|
| **CELLxGENE Census** (cellxgene.cziscience.com/census), accessed via the `cellxgene-census` Python SDK | Sole data source. |
| **He et al. 2022 fetal lung atlas**, dataset_id `3dc61ca1-ce40-46b6-8337-f27260fd9a03` (62,759 fetal-stage cells, 10 donors, 9–22 wpc, 10x 5′ v1) | Primary substrate; locked in [`notes/he2022_census_inventory.md`](notes/he2022_census_inventory.md). |
| `cellxgene-census` 1.17.0 + `scanpy` 1.11.5 + `anndata` 0.11.4 (+ transitive deps) | Whole stack. No other packages added. |

## Reproducing the v0 result

```sh
# Census inventory (substrate lock)
python3 scripts/tipcommit_census_inventory.py
# Q1 MVP (per-stage fractions + markers; 1–3 min runtime)
python3 scripts/tipcommit_q1_mvp.py
# Q1.1 refinement (state breakdown + curated markers; pure pandas)
python3 scripts/tipcommit_q1_refinement.py
```

Outputs land in `metadata/`. None of the scripts overwrite each
other's inputs.

## Repo layout

```
notes/
  evidence_map.md                project goal, Q1 scope, substrate
  status.md                      bootstrap state
  next_steps.md                  ordered 5-step plan
  he2022_census_inventory.md     substrate lock (Step 2 result)
  q1_mvp_design.md               MVP design contract (Step 3)
  q1_refinement.md               Q1.1 refinement note
  q1_summary.md                  v0 standing summary  ← start here
  q1_key_findings.md             8-bullet finding list  ← also here

scripts/
  tipcommit_census_inventory.py  Census substrate inventory
  tipcommit_q1_mvp.py            Q1 MVP (fractions + markers)
  tipcommit_q1_refinement.py     Q1.1 refinement (curation + restructure)

metadata/
  census_lung_datasets.csv         all lung-mentioning Census datasets
  census_fetal_lung_per_stage.csv  per-(dataset, stage) cell counts
  census_fetal_lung_celltypes.csv  cell-type distribution in chosen dataset
  q1_per_stage_fractions.csv       Q1 MVP wide format (3 thresholds)
  q1_long_format.csv               Q1 MVP figure-ready long format
  q1_markers_by_stage.csv          Q1 MVP top-20 airway/alveolar markers
  q1_state_breakdown_by_stage.csv  Q1.1 cleaner state breakdown
  q1_markers_by_stage_curated.csv  Q1.1 markers with gene_class flags
```

## Carry-over discipline from Gain

- Plan before code; small commits, push after each meaningful step.
- No generic framework if a project-specific script will do.
- No dependency unless it clearly reduces total complexity.
- Anti-overclaim: every note labels known biology vs inference vs
  hypothesis.
- No "future architecture" section longer than the actual working
  code.
