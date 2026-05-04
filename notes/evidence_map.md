# Evidence Map — TipCommit

## Project goal

Answer **Q1**: when and how does the bipotent SOX2 / SOX9
co-expressing tip cell commit to airway vs alveolar fate in **human**
fetal lung?

This project is the **expression-correlation route** chosen in
[Gain's next-project decision memo](https://github.com/DD-Ching/Gain/blob/main/notes/next_project_decision.md).
TipCommit is the successor to Gain's evidence-audit cycle, scoped
**strictly to Q1**. Q2 (continuous gradient vs discrete switch) is
explicitly deferred until Q1 has a working MVP.

## Scientific scope

The bipotent distal tip cell expresses **both SOX2 and SOX9** during
human fetal lung development. As lung morphogenesis proceeds, these
cells progressively commit to either:

- **proximal / airway lineage** (SOX2-positive, SOX9-negative;
  basal / secretory / ciliated cells), or
- **distal / alveolar lineage** (SOX9-positive, SOX2-negative;
  AT2 → AT1).

Q1 asks: **at what fetal stage** does this commitment happen, and
**what intermediate cell states** sit between the bipotent tip and
the committed daughters?

In scope:
- Human fetal lung scRNA-seq across multiple gestational stages
- SOX2 / SOX9 co-expression dynamics
- Marker gene shifts in committed-airway vs committed-alveolar
  daughter populations

Out of scope (for v0):
- Mouse data (different timing; cross-species comparison is v1+)
- Adult lung (Gain's audit covered this)
- Functional validation (lineage tracing, perturbation — wet lab)
- Q2 / proximal-distal axis continuity testing
- Other regulators (FGF10, WNT, BMP, SHH) beyond their role as
  marker genes
- Spatial transcriptomics analysis beyond using it as cross-check

## Known biology (textbook; cited in Gain)

- **NKX2-1** specifies lung epithelial identity from the embryonic
  endoderm.
- **SOX2** marks proximal airway progenitors; expressed in basal
  cells and persists in adult airway epithelium.
- **SOX9** marks distal tip / bud cells; canonical marker of the
  bipotent progenitor.
- **FGF10** (mesenchymal) signals to FGFR2b on distal epithelium,
  driving branching and maintaining the SOX9+ tip pool.
- **BMP, WNT, SHH** modulate the proximal-distal axis.

The textbook chain (stable):
```
NKX2-1 → SOX9+ distal tip → bipotent SOX2+SOX9+ → committed AT2 / AT1 (alveolar)
                                                → committed basal / secretory / ciliated (airway)
```

What is **not yet quantitatively settled in human**:
- The precise gestational stage at which commitment becomes
  irreversible.
- The cell-state continuum between bipotent and committed (single
  intermediate vs multiple).
- Whether SOX2 and SOX9 expression shift gradually or abruptly
  during commitment.

## Public resources

The expression substrate is mature; the audit substrate (Gain's
focus) is not used here.

### Primary substrate

- **CELLxGENE Census** (cellxgene.cziscience.com/census) —
  programmatic standardised access to human scRNA-seq via TileDB-
  SOMA. Python SDK `cellxgene-census` returns AnnData. Standardised
  cell-type, tissue, and developmental-stage ontology (CL, UBERON,
  HsapDv).

- **Human fetal lung atlas (He et al. 2022)** — *Nat Genet*. Staged
  human fetal lung scRNA-seq (5–22 weeks post-conception), with
  spatial transcriptomics. **The primary dataset** for Q1. Verify
  Census coverage at first inventory step.

### Secondary substrate (for cross-check)

- **Pan-fetal human cell atlas (Cao et al. 2020)** — *Science*. Lung
  subset gives proximal/distal progenitor coverage in a multi-organ
  context.

- **HLCA core (Sikkema et al. 2023)** — adult lung reference;
  useful as a landmark for "where do the committed cells end up?"
  but Q1 is a fetal question.

- **LungMAP human developmental scRNA-seq projects** — additional
  staged donors via the LungMAP Azul API. Use only if He 2022
  alone is insufficient.

### Tooling reused, not rebuilt

- **Scanpy** — `sc.pp` / `sc.tl` / `sc.pl`. Standard scRNA-seq
  workflow.
- **AnnData** — data container shared with Census + Scanpy.
- **scvi-tools** — out of scope for v0 (label transfer not needed
  for Q1 MVP).

## Initial deliverable (Q1 MVP)

The smallest useful Q1 result, to be scoped precisely in
`notes/next_steps.md` and locked in `notes/q1_mvp_design.md`:

- **Per fetal stage** (binned by gestational week or developmental
  stage label): fraction of epithelial cells in each of:
  - SOX2+ SOX9+ (bipotent)
  - SOX2+ SOX9− (airway-committed)
  - SOX2− SOX9+ (alveolar-committed)
  - SOX2− SOX9− (other / undefined)
- **Marker shifts**: at each stage, what genes distinguish committed-
  airway (SOX2+SOX9−) from committed-alveolar (SOX2−SOX9+) cells?
  Top 10–20 genes per direction.
- **Outputs**:
  - One CSV: per-stage fractions of the four cell classes.
  - One figure-ready table (long format): stage × class × count.
  - One marker-gene table: stage × top-N genes per direction.

What v0 does NOT do:
- Pseudotime / RNA-velocity / trajectory inference (defer to v1).
- Spatial cross-check (defer to v1+).
- Statistical inference on commitment-stage timing (descriptive
  only; confidence intervals are v1).
- Q2.

## Carry-over from Gain

Gain v0 (the predecessor evidence-audit project at
`github.com/DD-Ching/Gain`) is **frozen**. Its standing outputs
remain useful but TipCommit does not extend its audit cycle:

- **Substrate-gap diagnosis** — TipCommit accepts the public-data
  ChIP gap and pivots to scRNA-seq. The 5 NKX2-1 chromatin-
  supported standing distal candidates are preserved in Gain for
  future Hi-C / perturbation work.
- **Audit machinery** — Gain's 5-class evidence model is not
  applied here; expression analysis uses different idioms
  (Census query → AnnData → score → group-comparison).
- **Constraint discipline** — carries over verbatim:
  - No generic framework if a project-specific script will do.
  - No dependency unless it clearly reduces total complexity.
  - No "future architecture" section longer than the actual
    working code.
  - Plan before code.
  - Small commits, push after each meaningful step.
- **Anti-overclaim language** — carries over: distinguish known
  biology from inference from hypothesis; do not claim functional
  validation; do not blur correlation into causation.
