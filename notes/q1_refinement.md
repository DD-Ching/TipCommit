# Q1.1 Refinement Pass

**Date:** 2026-05-04
**Predecessor:** Q1 MVP — [`scripts/tipcommit_q1_mvp.py`](../scripts/tipcommit_q1_mvp.py),
committed at `a8ebd9d` along with the three `metadata/q1_*.csv` outputs.
**Script:** [`scripts/tipcommit_q1_refinement.py`](../scripts/tipcommit_q1_refinement.py)
**Outputs:**
- [`metadata/q1_state_breakdown_by_stage.csv`](../metadata/q1_state_breakdown_by_stage.csv)
- [`metadata/q1_markers_by_stage_curated.csv`](../metadata/q1_markers_by_stage_curated.csv)

This is a small refinement pass on top of the Q1 MVP. **No new
analysis surface, no Census re-fetch, no new dependencies.** The
script reads the existing per-stage and markers CSVs and restructures
them into a cleaner state breakdown and a flagged markers table.

## Goal A — clean state breakdown

[`q1_state_breakdown_by_stage.csv`](../metadata/q1_state_breakdown_by_stage.csv)
is a long-format table with one row per (stage, state). Six states:

| Family | State | Definition |
|---|---|---|
| candidate pool | `candidate_bipotent` | "epithelial cell of lung" with SOX2 > 0.5 **and** SOX9 > 0.5 |
| candidate pool | `transitioning_airway` | "epithelial cell of lung" with SOX2 > 0.5 only |
| candidate pool | `transitioning_alveolar` | "epithelial cell of lung" with SOX9 > 0.5 only |
| candidate pool | `undefined` | "epithelial cell of lung" with neither above 0.5 |
| annotated | `annotated_airway` | basal / club / multiciliated / secretory / respiratory tract / squamous |
| annotated | `annotated_alveolar` | AT1 / AT2 |

Each row carries `frac_of_total_epithelial` and (for the four
candidate-pool states) `frac_of_candidate_pool`. The candidate-pool
fraction is the within-pool view used for Goal C.

This table uses **only the moderate threshold** (>0.5). The full
threshold sweep stays in `q1_per_stage_fractions.csv` from the MVP.

## Goal B — curated markers

[`q1_markers_by_stage_curated.csv`](../metadata/q1_markers_by_stage_curated.csv)
keeps every row from `q1_markers_by_stage.csv` and adds two columns:

- `gene_class` — `biological` (default), `ribosomal` (matches
  `^(RPL|RPS|MRPL|MRPS)\d`), `mitochondrial` (matches `^MT-`),
  `housekeeping` (curated set: MALAT1, NEAT1, ZFAS1, XIST, NPM1,
  HNRNPA1, NACA, RACK1, GNAS, HSPA1A/B, HSP90AA1/AB1, HSBP1,
  MT2A/MT1X/MT1E/MT1G/MT1F, TMSB4X, TMSB10, ACTB, GAPDH, B2M,
  PPIA, HPRT1), or `unannotated` (Ensembl ID with no symbol).
- `recurrence_n_stages` — number of distinct stages this gene appears
  in for the same direction, counted only across the `biological`
  rows (so flagged genes don't inflate the count).

Class composition of the 160-row markers table:

| gene_class | rows | % |
|---|---:|---:|
| biological | 115 | 72% |
| ribosomal | 27 | 17% |
| housekeeping | 17 | 11% |
| unannotated | 1 | <1% |

**No rows are removed** — flagging only — so the original ranks are
still recoverable from this CSV alone.

### Repeated airway markers (biological, ≥2 of 4 successful stages)

| Gene | n_stages | best rank | stages | Note |
|---|---:|---:|---|---|
| AGR2 | 4 | 1 | 15,18,20,22 | secretory / anterior gradient |
| AGR3 | 4 | 1 | 15,18,20,22 | secretory / anterior gradient |
| CAPS | 4 | 3 | 15,18,20,22 | calcyphosine — multiciliated |
| HMGN3 | 3 | 3 | 15,18,22 | chromatin |
| S100A11 | 3 | 3 | 15,18,22 | calcium binding |
| IFT57 | 3 | 4 | 15,18,20 | intraflagellar transport |
| PLPP2 | 3 | 9 | 15,18,22 | lipid phosphatase |
| **SOX2** | **2** | **8** | **15,18** | **canonical airway TF** |
| ANXA1 | 2 | 4 | 15,18 | annexin |
| PERP | 2 | 5 | 15,18 | p53 effector / epithelial |
| FXYD3 | 2 | 7 | 15,18 | membrane modulator |
| EZR | 2 | 9 | 15,18 | apical cytoskeleton |
| DSTN | 2 | 3 | 15,18 | actin depolymerising |
| H2AJ | 2 | 13 | 15,18 | histone variant |
| TSPAN1 | 2 | 7 | 20,22 | tetraspanin |
| SMIM22 | 2 | 5 | 20,22 | small membrane integral |
| CIMAP1B | 2 | 5 | 20,22 | ciliary microtubule |
| CFAP144 | 2 | 10 | 20,22 | cilia / flagella associated |
| CFAP276 | 2 | 17 | 20,22 | cilia / flagella associated |
| CIMIP1 | 2 | 11 | 20,22 | cilia / microtubule |
| SPMIP6 | 2 | 11 | 20,22 | sperm/cilia microtubule |
| SPACA9 | 2 | 16 | 20,22 | sperm/cilia associated |
| DYNLL1 | 2 | 15 | 18,20 | dynein light chain |
| CALM1 | 2 | 8 | 15,22 | calmodulin |

Two clear sub-programs:
- **15–18 wpc airway shift**: secretory/anterior-gradient + apical/cytoskeletal
  (AGR2/3, ANXA1, EZR, DSTN, FXYD3, PERP) plus the airway TF **SOX2**.
- **20–22 wpc airway shift**: motile-cilia machinery dominates
  (CFAP144/276, CIMAP1B, CIMIP1, SPMIP6, SPACA9, TSPAN1, DYNLL1).

The two programs share AGR2/AGR3/CAPS at the top across all four stages.

### Repeated alveolar markers (biological, ≥2 of 4 successful stages)

| Gene | n_stages | best rank | stages | Note |
|---|---:|---:|---|---|
| **AGER** | 4 | 2 | 15,18,20,22 | **canonical AT1** |
| MYL9 | 4 | 1 | 15,18,20,22 | myosin light chain (also in AT1) |
| CLIC3 | 4 | 1 | 15,18,20,22 | chloride channel |
| SMARCA5 | 4 | 2 | 15,18,20,22 | chromatin remodeller |
| SPARC | 3 | 3 | 18,20,22 | matricellular |
| **CLDN18** | 2 | 7 | 15,18 | **alveolar tight junction** |
| **SFTPB** | 2 | 7 | 15,20 | **canonical AT2 surfactant** |
| BCAM | 2 | 8 | 15,18 | basal cell adhesion |

Plus stage-unique but biologically pointed hits:
- **NKX2-1** appears at 20 wpc (rank 5, alveolar direction) — distal-lung TF.
- CLDN6 (18), CPM (18), FOLR1 (18), MDK (18), AGRN (15), SDC1 (15),
  EMP2 (18), NREP (20), SERPINH1 (22), CRLF1 (15), RNASE1 (15),
  C19orf33 (15), S100A10 (18) all show up at single stages with
  alveolar / extracellular-matrix relevance.

Compared with the airway side, the alveolar side has **fewer
biological hits per stage** — most of the recurring rows are
dominated by ribosomal genes (RPL/RPS), which is why so few
non-flagged genes survive. The biological signal on the alveolar
side is real (AGER, SFTPB, CLDN18, NKX2-1) but lower-magnitude than
the airway side's secretory/ciliated programs.

## Goal C — trend robustness check

Same `candidate_bipotent` count (moderate threshold, >0.5 on both
SOX2 and SOX9) reported two ways:

| Stage | n | as % of total epithelial | as % of candidate_pool |
|---|---:|---:|---:|
| 9 wpc | 393 | 21.9% | 28.5% |
| 11 wpc | 499 | 24.4% | 35.4% |
| 15 wpc | 335 | 17.4% | 23.6% |
| 18 wpc | 86 | 4.8% | 6.9% |
| 20 wpc | 2 | 0.2% | 2.9% |
| 22 wpc | 1 | 0.1% | 1.5% |

**Verdict:** the decline is annotation-robust. Both denominators
trace the same shape — a small early bump 9 → 11 wpc, then a strong
and strictly monotonic decline 11 → 22 wpc. The drop is not an
artefact of the candidate pool itself shrinking relative to
committed cells (it does shrink, but the candidate-bipotent fraction
*within* that shrinking pool also collapses).

**Caveats stay the same:**
- The candidate pool is too small at 20–22 wpc (n = 70, 66) for the
  within-pool percentages to be statistically meaningful at those
  late stages (the n=2 and n=1 numerators sit at the edge of noise).
  The within-pool collapse is consistent with the trend but is
  documented at small n.
- 9–11 wpc still have **zero annotated_alveolar cells**, so the
  airway-vs-alveolar marker comparison cannot run there. The
  decline of `candidate_bipotent` itself is not affected (no
  marker DE is needed for it).
- `candidate_bipotent` remains a **candidate** state, not a proven
  lineage state. The expression-only labelling does not prove these
  cells are bipotent tip progenitors; it proves only that they
  co-express SOX2 and SOX9 above the moderate threshold inside the
  bulk-annotated "epithelial cell of lung" pool.

## What this refinement does **not** change

- The MVP's three CSVs are untouched. The refinement reads them and
  derives new shapes; it does not re-run normalisation, classification,
  or DE.
- The threshold sweep (relaxed/moderate/stringent) and the markers
  table both still live in their original MVP files for reference.
- No plotting, no Q2, no new packages added.
