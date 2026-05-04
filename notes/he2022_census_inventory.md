# He 2022 Census Inventory

**Run:** 2026-05-04 (UTC)
**Script:** [`scripts/tipcommit_census_inventory.py`](../scripts/tipcommit_census_inventory.py)
**Census version:** `stable` (the LTS snapshot at the time of run; the
`cellxgene_census` Python client picks the current LTS automatically)
**Outputs:**
- [`metadata/census_lung_datasets.csv`](../metadata/census_lung_datasets.csv) — all lung-mentioning datasets
- [`metadata/census_fetal_lung_per_stage.csv`](../metadata/census_fetal_lung_per_stage.csv) — per-(dataset, stage) cell counts
- [`metadata/census_fetal_lung_celltypes.csv`](../metadata/census_fetal_lung_celltypes.csv) — cell-type distribution in the chosen dataset

This note follows the rule from `notes/evidence_map.md`: separate
**confirmed**, **inferred**, and **uncertain** explicitly.

## Directly confirmed from CELLxGENE Census

These come straight from Census's `census_info.datasets` table and from
querying `census_data.homo_sapiens.obs`:

- **CELLxGENE Census contains a fetal lung scRNA-seq atlas matching
  the He 2022 specification.** The collection is titled
  *"A human fetal lung cell atlas uncovers proximal-distal gradients
  of differentiation and key regulators of epithelial fates"* — the
  paper's exact title.
- **The integrated "All cells" dataset is `3dc61ca1-ce40-46b6-8337-f27260fd9a03`.**
  62,759 cells in the fetal-stage subset (71,752 total cells reported
  by Census's per-dataset metadata; the difference reflects cells
  whose `tissue_general` filter or stage filter excluded them).
- **10 donors, single assay (10x 5′ v1), 100% primary data.**
- **Six fetal stage labels** (post-conception weeks; pseudoglandular
  through canalicular / saccular):

  | Stage | Cells |
  |---|---:|
  | 9th wpc | 8,280 |
  | 11th wpc | 9,531 |
  | 15th wpc | 12,071 |
  | 18th wpc | 13,170 |
  | 20th wpc | 7,490 |
  | 22nd wpc | 12,217 |
  | **Total** | **62,759** |

- **75 distinct cell types observed** in this dataset.
- **9,517 of 62,759 cells (15%) match epithelial keywords** (epithel /
  AT1 / AT2 / alveol / basal / secret / ciliat / club / goblet /
  bud tip / tip cell). Per-type breakdown:

  | Cell type | Cells | Q1 role |
  |---|---:|---|
  | epithelial cell of lung | 5,595 | likely bipotent / undifferentiated tip-like — the SOX2/SOX9 substrate |
  | lung multiciliated epithelial cell | 1,641 | committed airway (ciliated) |
  | respiratory tract epithelial cell | 991 | airway broadly |
  | pulmonary alveolar type 2 cell | 382 | committed alveolar (AT2) |
  | pulmonary alveolar type 1 cell | 283 | committed alveolar (AT1) |
  | basal cell | 278 | committed airway (basal) |
  | club cell | 208 | committed airway (secretory / club) |
  | lung secretory cell | 128 | committed airway (secretory) |
  | squamous epithelial cell | 11 | likely metaplastic / outlier |

- **Across the broader Census slice (`tissue_general == 'lung'`):**
  - 10,181,125 lung cells total (across human + adult + fetal)
  - 92 distinct datasets
  - 1,401 distinct donors
  - 194 distinct development_stage labels (most are adult-year-old
    bins; 17 datasets contain fetal/embryonic/wpc-staged cells).

- **Other fetal-lung-mentioning datasets exist but are NOT the right
  substrate for Q1:**
  - `350237e0-9f48-4cbd-9140-3b44495549f3` ("Fetal lung + Pan-fetal
    immune"): 544k cells but mostly immune (only 869 epithelial-
    keyword matches). Wrong substrate for SOX2/SOX9 epithelial
    progenitor work.
  - Per-pcw subsets (`b3...`, `c9...`, `2a...`, `30...`, `47...`):
    individual donor / sample slices, ~5,000 cells each. The "All
    cells" integrated dataset (3dc61ca1) is the integration of these
    plus other donors. Use the integrated set, not the per-donor
    splits.
  - Organoid dataset (`4023a2bc...`): in the same collection but
    organoid-derived, not primary tissue. Useful as a v1+
    cross-check.

## Inferred (not directly verified at this step)

- **The "epithelial cell of lung" category (5,595 cells, 59% of
  epithelial) is most likely the bipotent / undifferentiated /
  tip-progenitor population.** It is the catch-all for epithelial
  cells that the upstream paper's annotation did not classify as
  AT1, AT2, basal, club, ciliated, etc. SOX2 and SOX9 expression
  scoring should pull the bipotent vs committed split out of this
  category. **To be verified at MVP-design step** by inspecting
  SOX2/SOX9 expression distributions within "epithelial cell of
  lung".
- **The 11–22 wpc range covers the pseudoglandular → canalicular
  → early saccular transition.** Q1's bipotent-to-committed window
  should be visible across these stages. **Inference, not
  Census-confirmed**: pseudoglandular stage is roughly 5–17 wpc and
  canalicular 16–26 wpc per textbook definitions; the dataset's
  stages span both.
- **The single-assay 10x 5′ v1 platform makes this dataset
  internally consistent** for differential expression — no
  cross-assay batch correction needed.

## Uncertain / to-be-resolved at the MVP-design step

- **SOX2 and SOX9 expression thresholds** for defining "+/−" cells.
  Will require inspecting the actual expression distribution in the
  epithelial subset (likely `> 0` in log-normalised space, but
  potentially a percentile-based cutoff if there is a long zero-
  inflated tail).
- **Whether the per-pcw donor distribution is balanced enough** to
  trust per-stage fractions. With 10 donors total across 6 stages,
  some stages have ~1–2 donors. Per-stage estimates are descriptive,
  not statistical.
- **Whether "epithelial cell of lung" is a single homogeneous
  bipotent pool or contains substructure** (e.g., proximal vs distal
  tip, cycling vs quiescent). Sub-clustering is v1 work; v0 treats
  it as one pool.
- **Whether the upstream paper's authors are "He et al. 2022" exactly
  or another author group.** The collection title matches papers in
  the proximal-distal gradient series; the exact citation to use is
  retrievable from the Census collection metadata (`gh` /
  cellxgene.cziscience.com lookup) or the dataset's `dataset_h5ad_path`
  publication metadata. Practical impact for v0 is zero — the
  substrate is correct regardless of which exact 2022/2024 atlas
  paper deposited it.

## Decision: can the Q1 MVP proceed on the He 2022 dataset directly?

**Yes.** Dataset `3dc61ca1-ce40-46b6-8337-f27260fd9a03` provides:

- ✅ Sufficient cell count (62,759 with 9,517 epithelial)
- ✅ Multi-stage coverage (6 fetal stages, 9–22 wpc — the relevant
  developmental window)
- ✅ Both bipotent-like (epithelial cell of lung, 5,595) and
  committed (alveolar AT1+AT2 = 665; airway basal+club+secretory+
  multiciliated = 2,255) populations annotated
- ✅ Single assay (10x 5′ v1) — no cross-platform batch issues
- ✅ Programmatic access via `cellxgene_census` Python SDK with
  standardised cell-type and stage ontologies
- ✅ 100% primary data (no duplicates / cross-dataset overlap)

No fallback needed. The MVP design step (`notes/q1_mvp_design.md`,
the next deliverable per `notes/next_steps.md`) can lock the Census
slice as:

```python
value_filter = (
    "tissue_general == 'lung' and "
    "dataset_id == '3dc61ca1-ce40-46b6-8337-f27260fd9a03'"
)
```

with cell-type filter applied downstream in pandas / scanpy.

## Reproducing this inventory

```sh
python3 scripts/tipcommit_census_inventory.py
```

Requires the `cellxgene-census`, `scanpy`, `anndata` stack
(installed at the start of step 2; documented in
[`notes/status.md`](status.md)). Runtime ~ 60–120 s (the
`obs.read` of the human-lung slice is the main cost). No new HTTP
beyond Census (which the SDK handles internally).
