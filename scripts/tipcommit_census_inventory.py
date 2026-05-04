#!/usr/bin/env python3
"""tipcommit census-inventory - locate the He et al. 2022 fetal lung dataset
(or its closest analogue) in CELLxGENE Census and report substrate coverage.

Step 2 of TipCommit's plan (notes/next_steps.md). Stops at inventory --
does NOT fetch expression data, does NOT begin the Q1 MVP.

Outputs:
- console summary
- metadata/census_lung_datasets.csv      (all lung-mentioning datasets)
- metadata/census_fetal_lung_per_stage.csv  (cell counts per dataset x stage)
- metadata/census_fetal_lung_celltypes.csv  (cell types in the chosen dataset)
"""
from __future__ import annotations

import sys
from pathlib import Path

import cellxgene_census  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "metadata"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Search terms used to identify candidate fetal lung datasets in titles
LUNG_KEYWORDS = ["lung"]
FETAL_KEYWORDS = ["fetal", "fetus", "embryonic", "developmental", "development"]
HE_KEYWORDS = ["he,", "he et", "he2022", "branching morphogenesis"]


def section(title: str) -> None:
    print(f"\n=== {title} ===")


def main() -> int:
    section("Opening Census (LTS 'stable')")
    with cellxgene_census.open_soma(census_version="stable") as census:
        # ---- Step A: lung-mentioning datasets ----
        section("Lung-mentioning datasets")
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
        print(f"Census datasets total: {len(datasets):,}")

        title_lower = datasets["dataset_title"].fillna("").str.lower()
        coll_lower = datasets["collection_name"].fillna("").str.lower()
        lung_mask = (
            title_lower.str.contains("|".join(LUNG_KEYWORDS), regex=True)
            | coll_lower.str.contains("|".join(LUNG_KEYWORDS), regex=True)
        )
        lung_datasets = datasets[lung_mask].copy()
        print(f"Lung-mentioning datasets: {len(lung_datasets)}")

        lung_datasets.to_csv(OUT_DIR / "census_lung_datasets.csv", index=False)
        print(f"  wrote metadata/census_lung_datasets.csv")

        # Look for fetal/developmental subset
        title_lower_lung = lung_datasets["dataset_title"].fillna("").str.lower()
        coll_lower_lung = lung_datasets["collection_name"].fillna("").str.lower()
        fetal_mask = (
            title_lower_lung.str.contains("|".join(FETAL_KEYWORDS), regex=True)
            | coll_lower_lung.str.contains("|".join(FETAL_KEYWORDS), regex=True)
        )
        fetal_lung_datasets = lung_datasets[fetal_mask].copy()
        print(f"\nFetal/developmental lung-mentioning datasets: {len(fetal_lung_datasets)}")
        for _, row in fetal_lung_datasets.iterrows():
            print(f"  dataset_id: {row['dataset_id']}")
            print(f"    title:      {row['dataset_title'][:120]}")
            print(f"    collection: {row['collection_name'][:120]}")
            print(f"    cells:      {row.get('dataset_total_cell_count', '?')}")
            print()

        # Look for He 2022 specifically
        section("Looking for He et al. 2022 by keyword")
        he_mask = title_lower_lung.str.contains("|".join(HE_KEYWORDS), regex=True)
        he_hits = lung_datasets[he_mask]
        if len(he_hits) > 0:
            print(f"He-keyword hits: {len(he_hits)}")
            for _, row in he_hits.iterrows():
                print(f"  {row['dataset_id']}  {row['dataset_title'][:120]}")
        else:
            print("No direct title match for He 2022 keywords. ")
            print("Will rely on the 'fetal lung' subset above as candidates.")

        # ---- Step B: query human lung obs for substrate coverage ----
        section("Querying obs for human lung tissue (this can take ~30-60 s)")
        obs_lung = census["census_data"]["homo_sapiens"].obs.read(
            value_filter="tissue_general == 'lung'",
            column_names=[
                "dataset_id", "development_stage", "cell_type", "tissue",
                "donor_id", "is_primary_data", "assay",
            ],
        ).concat().to_pandas()
        print(f"Total human lung cells in Census: {len(obs_lung):,}")
        print(f"Distinct datasets:                 {obs_lung['dataset_id'].nunique()}")
        print(f"Distinct donors:                   {obs_lung['donor_id'].nunique()}")

        section("Development-stage breakdown across all human lung")
        stages_all = obs_lung["development_stage"].value_counts()
        stages_all_nonzero = stages_all[stages_all > 0]
        print(f"Distinct stages: {len(stages_all_nonzero)}")
        print("Top 30 stages by cell count:")
        print(stages_all_nonzero.head(30).to_string())

        # ---- Step C: filter to fetal-stage subset ----
        section("Fetal-stage subset")
        fetal_kw_pattern = "|".join(FETAL_KEYWORDS + ["week", "wpc"])
        # development_stage is a pandas Categorical; cast to str before .str.contains
        stage_str = obs_lung["development_stage"].astype(str).str.lower()
        fetal_stage_mask = stage_str.str.contains(
            fetal_kw_pattern, regex=True, na=False,
        )
        obs_fetal = obs_lung[fetal_stage_mask]
        print(f"Fetal-stage human lung cells: {len(obs_fetal):,}")
        print(f"Distinct datasets:             {obs_fetal['dataset_id'].nunique()}")
        print(f"Distinct donors:               {obs_fetal['donor_id'].nunique()}")

        # Per-dataset summary in fetal subset
        section("Per-dataset summary (fetal stages only)")
        by_dataset = obs_fetal.groupby("dataset_id").agg(
            n_cells=("cell_type", "size"),
            n_donors=("donor_id", "nunique"),
            n_celltypes=("cell_type", "nunique"),
            n_stages=("development_stage", "nunique"),
        ).reset_index().sort_values("n_cells", ascending=False)
        print(by_dataset.head(20).to_string(index=False))

        # Save per-stage breakdown (drop zero-count combinations from
        # the Categorical cross-product to keep the CSV small)
        if len(obs_fetal) > 0:
            per_stage = (
                obs_fetal.groupby(["dataset_id", "development_stage"], observed=True)
                .size()
                .reset_index(name="n_cells")
            )
            per_stage = per_stage[per_stage["n_cells"] > 0]
            per_stage = per_stage.sort_values(
                ["dataset_id", "n_cells"], ascending=[True, False],
            )
            per_stage.to_csv(OUT_DIR / "census_fetal_lung_per_stage.csv", index=False)
            print(f"\nwrote metadata/census_fetal_lung_per_stage.csv "
                  f"({len(per_stage)} nonzero rows)")

        # ---- Step D: deep dive on chosen He 2022 atlas dataset ----
        # The script's first run revealed two distinct candidate datasets:
        #   350237e0...  "Fetal lung + Pan-fetal immune" — 544k cells but
        #                mostly immune (not what Q1 needs)
        #   3dc61ca1...  "All cells" from the proximal-distal-gradient atlas
        #                — 62k cells, 75 cell types, 10 donors, 6 stages.
        # The He 2022 atlas is the right substrate for SOX2/SOX9 epithelial
        # progenitor work; the immune dataset is not.
        HE_ATLAS_DATASET_ID = "3dc61ca1-ce40-46b6-8337-f27260fd9a03"
        if len(obs_fetal) == 0:
            print("\nNo fetal-stage human lung cells found in Census.")
            print("This is a meaningful negative for TipCommit's plan.")
            return 0

        if HE_ATLAS_DATASET_ID in set(obs_fetal["dataset_id"]):
            chosen_did = HE_ATLAS_DATASET_ID
            print(f"\nUsing He-2022-style atlas dataset {chosen_did} for deep dive "
                  f"(preferred over the largest immune-cell-only dataset).")
        else:
            chosen_did = obs_fetal["dataset_id"].value_counts().idxmax()
            print(f"\nHe atlas dataset {HE_ATLAS_DATASET_ID} not found in fetal "
                  f"subset; falling back to largest: {chosen_did}")
        section(f"Deep dive on chosen fetal lung dataset: {chosen_did}")
        sub = obs_fetal[obs_fetal["dataset_id"] == chosen_did]
        print(f"  cells:   {len(sub):,}")
        print(f"  donors:  {sub['donor_id'].nunique()}")
        print(f"  assays:  {sub['assay'].dropna().unique().tolist()[:5]}")
        print(f"  primary-data fraction: {sub['is_primary_data'].mean():.2f}")

        ds_meta = datasets[datasets["dataset_id"] == chosen_did]
        if len(ds_meta) > 0:
            r = ds_meta.iloc[0]
            print(f"\n  Title:      {r['dataset_title']}")
            print(f"  Collection: {r['collection_name']}")
            print(f"  Total cells in dataset: {r.get('dataset_total_cell_count', '?')}")

        print(f"\n  Stage distribution (observed only):")
        stage_counts = sub["development_stage"].value_counts()
        print(stage_counts[stage_counts > 0].to_string())

        ct_counts = sub["cell_type"].value_counts()
        ct_counts_nonzero = ct_counts[ct_counts > 0]
        print(f"\n  Cell-type distribution (top 30 of "
              f"{len(ct_counts_nonzero)} observed):")
        print(ct_counts_nonzero.head(30).to_string())

        # Save observed cell-type list for the chosen dataset (drop the
        # Categorical-zeroes from other tissues' cell types).
        ct_counts_nonzero.to_csv(OUT_DIR / "census_fetal_lung_celltypes.csv",
                                 header=["n_cells"])
        print(f"\nwrote metadata/census_fetal_lung_celltypes.csv "
              f"({len(ct_counts_nonzero)} observed cell types)")

        # ---- Step E: targeted check for SOX2 / SOX9 epithelial coverage ----
        section("Epithelial cell type coverage in chosen dataset")
        epi_keywords = ["epithel", "AT1", "AT2", "alveol", "bud tip", "tip cell",
                       "basal", "secret", "ciliat", "club", "goblet"]
        epi_pattern = "|".join(epi_keywords)
        ct_str = sub["cell_type"].astype(str).str.lower()
        epi_mask = ct_str.str.contains(epi_pattern, regex=True, na=False)
        epi_cells = sub[epi_mask]
        print(f"Epithelial-keyword cells: {len(epi_cells):,} of {len(sub):,} total")
        print(f"Epithelial cell types in this dataset:")
        epi_counts = epi_cells["cell_type"].value_counts()
        print(epi_counts[epi_counts > 0].to_string())

    print("\nInventory complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
