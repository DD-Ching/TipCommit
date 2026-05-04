# Status — TipCommit

**Date:** 2026-05-04
**Phase:** bootstrap (planning notes)
**Predecessor:** [Gain v0.1.0](https://github.com/DD-Ching/Gain/releases/tag/v0.1.0) (frozen)

## Starting state

- Repo started **fresh** at `/Users/ddh/Downloads/TipCommit`.
  Separate from Gain per the
  [next-project decision](https://github.com/DD-Ching/Gain/blob/main/notes/next_project_decision.md).
- **No code yet — planning phase.** Three planning notes
  (`notes/evidence_map.md`, `notes/status.md`, `notes/next_steps.md`)
  before any implementation.
- **Strictly Q1 only.** Q2 (proximal-distal axis continuity) is
  deferred.
- **Gain's standing outputs are preserved at the predecessor repo**
  but not extended here. TipCommit can cite Gain; it does not
  re-derive Gain's findings.

## Tooling check

- `git` 2.50.1 — available
- `gh` 2.89.1 — authenticated as `DD-Ching`
- `python3` 3.10.12 — available
- **Python single-cell stack: INSTALLED 2026-05-04.**
  - `cellxgene-census` 1.17.0
  - `scanpy` 1.11.5
  - `anndata` 0.11.4
  - Plus transitive deps (numpy, pandas, scipy, tiledbsoma, pyarrow).
  - Install command used: `python3 -m pip install cellxgene-census scanpy anndata`
    (into the existing miniforge base env; no venv).

This is a **deliberate dep-add** vs Gain's stdlib-only constraint —
acknowledged in `notes/next_project_decision.md`. The new scope
discipline rule for TipCommit is "no packages beyond the
cellxgene-census + scanpy + anndata stack" (and their transitive
deps, e.g., numpy, pandas, scipy).

## What exists in the repo right now

```
notes/
  evidence_map.md           project goal, Q1 scope, substrate, MVP outline
  status.md                 this file
  next_steps.md             ordered actions with exit conditions
  he2022_census_inventory.md  step-2 inventory result (substrate locked)
metadata/
  census_lung_datasets.csv         all lung-mentioning Census datasets
  census_fetal_lung_per_stage.csv  per-(dataset, stage) cell counts
  census_fetal_lung_celltypes.csv  cell-type distribution in chosen dataset
scripts/
  tipcommit_census_inventory.py    Step 2 inventory script
```

Step 2 (Census inventory) is **complete**. Step 3 (Q1 MVP design)
and Step 4 (MVP implementation) are the next deliverables, per
`notes/next_steps.md`.

## What is decided

- **Project goal:** Q1 only (bipotent SOX2/SOX9 commitment timing
  in human fetal lung).
- **Substrate:** CELLxGENE Census, primarily the He et al. 2022
  fetal lung dataset.
- **Method idiom:** Census query → AnnData → expression scoring →
  per-stage fractions + marker comparison. No pseudotime, no
  trajectory inference, no scvi-tools in v0.
- **Output shape:** one CSV (per-stage fractions) + one figure-
  ready table (long format) + one marker-gene table.
- **Repo name:** `TipCommit` (descriptive — distal tip cell
  commitment).

## What is *not* decided yet

- **GitHub remote:** repository on github.com/DD-Ching/TipCommit will
  be created when this scaffold is committed. Until then, local-only.
- **Exact stage binning:** depends on what gestational-stage labels
  He 2022 / Census use. Locked at the inventory step in
  `next_steps.md`.
- **SOX2+/SOX9+ thresholds:** depends on observed expression
  distributions. Locked at the MVP-design step.
- **Differential-expression method:** Scanpy's
  `sc.tl.rank_genes_groups` defaults are the starting point; tune
  if results are noisy.

## Carry-over rules from Gain v0

- Plan before code.
- Small commits, push after each meaningful step.
- No generic framework if a project-specific script will do.
- No "future architecture" section longer than the actual working
  code.
- Anti-overclaim: separate **known biology** from **inference** from
  **hypothesis** in every note.
- The sciVI-tools-equivalent tool (any heavy single-cell library) is
  added only when a concrete script demands it.
