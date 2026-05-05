# TipCommit

> **Bipotent SOX2 / SOX9 commitment in human fetal lung
> epithelium — Q1, Q1.5b cross-dataset replication, and Q2 axis-
> continuity v0 cycles complete.** Census + Scanpy. Stdlib + the
> Census/Scanpy/AnnData stack only.

A small expression-correlation project answering both halves of a
single biological question on the He et al. 2022 fetal lung atlas
(via CELLxGENE Census), with a directional cross-dataset check on
Cao et al. 2020:

> **Q1 — when:** when do bipotent SOX2 / SOX9 co-expressing tip
> cells commit to airway vs alveolar fate in human fetal lung?
>
> **Q2 — how:** is the bipotent → committed transition continuous
> (gradient) or discrete (switch)?

TipCommit is the successor to
[Gain v0.1.0](https://github.com/DD-Ching/Gain/releases/tag/v0.1.0)
(an evidence-audit / substrate-gap repo whose chain stopped at the
public-data ChIP gap). TipCommit picks up the
expression-correlation route documented in
[Gain's next-project decision memo](https://github.com/DD-Ching/Gain/blob/main/notes/next_project_decision.md).

## Project status

**Q1 v0, Q1.5b cross-dataset replication, and Q2 v0 cycles all
complete (2026-05-05).** Standing outputs frozen in this repo
phase. No further passes planned without a new explicit question.

### Start here (in order)

For the **Q1 result** ("when does commitment happen?"):

1. **[`notes/q1_key_findings.md`](notes/q1_key_findings.md)** —
   eight bullets, two-minute read.
2. **[`notes/q1_summary.md`](notes/q1_summary.md)** — full v0
   narrative.
3. **[`notes/q1_refinement.md`](notes/q1_refinement.md)** — Q1.1
   refinement (cleaner state breakdown + curated markers).
4. **[`notes/q15b_replication.md`](notes/q15b_replication.md)** —
   Q1.5b cross-dataset replication on Cao 2020.

For the **Q2 result** ("is the transition gradient or switch?"):

5. **[`notes/q2_key_findings.md`](notes/q2_key_findings.md)** —
   eight bullets, two-minute read.
6. **[`notes/q2_summary.md`](notes/q2_summary.md)** — full v0
   narrative including the Q2.1 ambiguity-resolution pass.
7. **[`notes/q21_candidate_pool_only.md`](notes/q21_candidate_pool_only.md)** —
   the within-pool follow-up that overturned the late-stage
   switch reading.

Design contracts and substrate inventories sit alongside their
result notes — see *Repo layout* below.

## What this repo establishes

### Q1 — when

1. **A robust per-stage trend in human fetal lung.** The candidate
   SOX2 / SOX9 co-expressing fraction inside the He 2022 atlas's
   candidate epithelial pool collapses from ~24% at 11 wpc to ~0.1%
   by 22 wpc at the moderate threshold. The trajectory is the same
   under three independent thresholds and under both denominator
   choices (total epithelial vs candidate pool only).
2. **Coherent stage-structured marker programs.** Airway-shifted
   markers organise into a 15–18 wpc secretory / apical program
   (AGR2, AGR3, CAPS, ANXA1, EZR, PERP) and a 20–22 wpc
   motile-cilia program (CFAP144/276, CIMAP1B, CIMIP1, TSPAN1).
   Alveolar-shifted markers recover the canonical AT1 / AT2 set
   (AGER, SFTPB, CLDN18) plus NKX2-1 at 20 wpc.

### Q1.5b — replication

3. **Directional support from an independent dataset.** The same
   pipeline on Cao et al. 2020 (Science 2020, sci-RNA-seq3,
   different lab, different annotation pipeline) shows: same
   direction (decline 12 → 17 wpc), much smaller magnitude
   (~50–100× lower fractions, attributable to documented assay
   sensitivity gap), substituted but biologically equivalent
   markers (SCGB3A2 airway; SFTPC + ETV5 + SLC34A2 alveolar). Cao
   covers only 12–17 wpc and so cannot test He's early bump or
   late collapse.

### Q2 — how

4. **Within the candidate pool, the bipotent → committed
   transition is gradient-like at every observable stage.** The
   within-pool bimodality coefficient on a per-cell commitment
   score (multi-gene airway − alveolar signature) stays at 0.31–
   0.43 across 9 → 22 wpc — uniformly below the 0.555 bimodality
   threshold. Co-commitment (cells expressing both programs at
   once) peaks at 34% inside the pool at 15 wpc — a textbook
   gradient signature.
5. **The all-epithelial-level mixed pattern is a composition
   effect, not within-pool hardening.** When AT1 / AT2 / basal /
   club / multiciliated / secretory cells are included, BC flips
   to switch-like (>0.7) at 15+ wpc. When restricted to the
   candidate pool only, the late-stage switch disappears (15 wpc
   BC drops from 0.71 to 0.38). The all-cells mixed pattern is the
   joint result of (Q1) pool depletion plus (Q2) a steady within-
   pool gradient — not anyone hardening within the pool.

### Method discipline

6. **A reproducible MVP + refinement pipeline.** Five scripts
   running end-to-end against a public Census release on the
   `cellxgene-census` + `scanpy` + `anndata` stack with no other
   dependencies. No plotting infrastructure; figure-ready CSVs
   only.

## What this repo is *not*

- **Not a lineage-tracing claim.** SOX2 + SOX9 co-expression and
  the multi-gene signatures define candidate / expression-level
  states — not proven lineage states. `candidate_bipotent` stays
  a candidate label throughout.
- **Not a causal or regulatory claim.** Q1 + Q2 do not test what
  drives the choice between fates; the Gain-era public-data ChIP
  gap still applies.
- **Not a quantitative cross-dataset replication.** Cao 2020
  supports direction but not magnitude (Q1.5b documented its
  assay-sensitivity limit). Q2 used Cao directionally only.
- **Not a perturbation, spatial, or pseudotime analysis.** The
  ordering is the dataset's `development_stage` label; no
  trajectory inference.

## What remains open

The substrate-bound questions inside the bipotent → committed
window are answered to the level a public scRNA-seq atlas can
reach. What is genuinely *not* answered:

- **What drives the choice.** Whether SOX2, SOX9, FGF, WNT, BMP,
  or another input regulator causes a given cell to commit
  airway vs alveolar — a perturbation question, out of substrate
  scope for this repo.
- **Lineage proof of any state.** Whether `candidate_bipotent`
  cells are functionally bipotent — needs lineage tracing or
  live imaging.
- **Sub-structure inside the candidate pool (Q1.5a — deferred).**
  Whether `epithelial cell of lung` contains proximal-vs-distal-tip
  or cycling / quiescent sub-populations.
- **Cross-species comparison.** Whether the same gradient-in-pool
  pattern holds in mouse fetal lung at the equivalent stages.
- **Donor-level variance.** With 10 He donors across 6 stages and
  11 Cao donors across 5, some stages are 1–2 donors. v0 reports
  fractions and shape diagnostics descriptively only.
- **Spatial cross-check.** Whether the gradient-like in-pool
  distribution is also spatial-gradient-like along the proximal-
  distal axis (would need spatial transcriptomics).

## Recommended next direction

The remaining open items above all require either (a) a different
data modality (perturbation, lineage, spatial) or (b) a non-trivial
new analysis question. None is a natural one-pass extension of the
existing TipCommit pipeline. Plausible v1+ moves:

- **Mouse cross-species replication** (Tabula Muris Senis fetal
  subsets, Han 2020 mouse cell atlas) — same pipeline, new
  substrate, tests whether the within-pool gradient pattern is
  evolutionarily conserved.
- **Multi-signature ensemble** — re-run Q2 with several alternative
  airway / alveolar / distal-tip composites and report the
  bimodality reading as a distribution, not a single number.
- **Q3 — proximal-distal axis spatial check** using a fetal-lung
  spatial transcriptomics dataset; tests whether the in-pool
  gradient is also a spatial gradient.

Each is a separate-design proposition. None is in scope for this
repo's current phase.

## Substrate

| Resource | Use |
|---|---|
| **CELLxGENE Census** (cellxgene.cziscience.com/census, version `2025-11-08`), accessed via the `cellxgene-census` Python SDK | Sole data source. |
| **He et al. 2022 fetal lung atlas**, dataset_id `3dc61ca1-ce40-46b6-8337-f27260fd9a03` (62,759 fetal-stage cells, 10 donors, 9–22 wpc, 10x 5′ v1) | **Primary substrate** for Q1, Q2; locked in [`notes/he2022_census_inventory.md`](notes/he2022_census_inventory.md). |
| **Cao et al. 2020 fetal cell atlas (1M subset)**, dataset_id `fa27492b-82ff-4ab7-ac61-0e2b184eee67` (53,429 lung cells, 11 donors, 12–17 wpc, sci-RNA-seq3) | **Secondary substrate** for Q1.5b replication and Q2 directional support; locked in [`notes/q15b_replication_design.md`](notes/q15b_replication_design.md). |
| `cellxgene-census` 1.17.0 + `scanpy` 1.11.5 + `anndata` 0.11.4 (+ transitive deps: numpy, pandas, scipy, tiledbsoma) | Whole stack. **No other packages added across the entire repo.** |

## Reproducing the v0 result

```sh
# Substrate inventory (Q1 prerequisite; locks He 2022)
python3 scripts/tipcommit_census_inventory.py
# Q1 MVP (per-stage fractions + markers; ~1-3 min)
python3 scripts/tipcommit_q1_mvp.py
# Q1.1 refinement (state breakdown + curated markers; pure pandas)
python3 scripts/tipcommit_q1_refinement.py
# Q1.5b replication (Cao 2020 direction check; ~2-4 min)
python3 scripts/tipcommit_q15b_replication.py
# Q2 MVP (multi-gene scoring + bimodality diagnostics; ~10-25 min;
#         dominated by Census fetch of two large slices)
python3 scripts/tipcommit_q2_mvp.py
# Q2.1 ambiguity-resolution pass (pure pandas restructure)
python3 scripts/tipcommit_q21_candidate_pool_only.py
```

Outputs land in `metadata/`. None of the scripts overwrite each
other's inputs. The Q2.1 script depends on Q2 MVP's per-cell
scores CSV; everything else can run independently.

## Repo layout

```
notes/
  evidence_map.md                 project goal, Q1 + Q2 scope, substrate
  status.md                       bootstrap state
  next_steps.md                   ordered 5-step plan
  he2022_census_inventory.md      Q1 substrate lock
  q1_mvp_design.md                Q1 MVP design contract
  q1_refinement.md                Q1.1 refinement note
  q1_summary.md                   Q1 v0 standing summary
  q1_key_findings.md              Q1 8-bullet findings
  q15b_replication_design.md      Q1.5b design contract
  q15b_replication.md             Q1.5b result writeup
  q2_design.md                    Q2 design contract
  q2_mvp.md                       Q2 MVP result writeup
  q21_candidate_pool_only.md      Q2.1 ambiguity-resolution pass
  q2_summary.md                   Q2 v0 standing summary
  q2_key_findings.md              Q2 8-bullet findings

scripts/
  tipcommit_census_inventory.py   Census substrate inventory
  tipcommit_q1_mvp.py             Q1 MVP (fractions + markers)
  tipcommit_q1_refinement.py      Q1.1 refinement (curation + restructure)
  tipcommit_q15b_replication.py   Q1.5b Cao 2020 replication
  tipcommit_q2_mvp.py             Q2 MVP (multi-gene scoring + diagnostics)
  tipcommit_q21_candidate_pool_only.py  Q2.1 within-pool diagnostic pass

metadata/
  census_lung_datasets.csv         all lung-mentioning Census datasets
  census_fetal_lung_per_stage.csv  per-(dataset, stage) cell counts
  census_fetal_lung_celltypes.csv  cell-type distribution in chosen dataset
  q1_per_stage_fractions.csv       Q1 MVP wide format (3 thresholds)
  q1_long_format.csv               Q1 MVP figure-ready long format
  q1_markers_by_stage.csv          Q1 MVP top-20 airway/alveolar markers
  q1_state_breakdown_by_stage.csv  Q1.1 cleaner state breakdown
  q1_markers_by_stage_curated.csv  Q1.1 markers with gene_class flags
  q15b_replication_per_stage_fractions.csv  Cao 2020 wide-format replication
  q15b_replication_long_format.csv          Cao 2020 long format with denom
  q15b_replication_markers_by_stage.csv     Cao 2020 markers (curated)
  q2_per_cell_scores.csv           Q2 per-cell signature scores (both substrates)
  q2_per_stage_diagnostics.csv     Q2 He 2022 per-stage BC + gap + co-commit
  q2_cao_secondary_diagnostics.csv Q2 Cao 2020 direction-only diagnostics
  q2_joint_density.csv             Q2 5x5 joint-density grids per stage
  q21_candidate_pool_diagnostics.csv  Q2.1 within-pool diagnostics
```

## Carry-over discipline from Gain

- Plan before code; small commits, push after each meaningful step.
- No generic framework if a project-specific script will do.
- No dependency unless it clearly reduces total complexity.
- Anti-overclaim: every note labels known biology vs inference vs
  hypothesis.
- No "future architecture" section longer than the actual working
  code.
