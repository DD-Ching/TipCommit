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
- **Python single-cell stack: NOT YET INSTALLED.**
  - `scanpy` — required for v0 MVP
  - `cellxgene-census` — required for Census query
  - `anndata` — required for in-memory representation
  - Install pending; will be done at the start of MVP
    implementation, not earlier (per "no dependency unless it
    clearly reduces total complexity" — install when needed).

This is a **deliberate dep-add** vs Gain's stdlib-only constraint —
acknowledged in `notes/next_project_decision.md`. The new scope
discipline rule for TipCommit is "no packages beyond the
cellxgene-census + scanpy + anndata stack" (and their transitive
deps, e.g., numpy, pandas, scipy).

## What exists in the repo right now

```
notes/
  evidence_map.md     project goal, Q1 scope, substrate, MVP outline
  status.md           this file
  next_steps.md       ordered actions with exit conditions
metadata/             (empty; will receive Census query results, MVP CSVs)
scripts/              (empty; will receive Q1 MVP script)
```

That's it. No code, no dependency manifest, no CI. Those land when
there is a concrete reason for them.

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
