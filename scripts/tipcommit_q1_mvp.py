#!/usr/bin/env python3
"""tipcommit q1-mvp - first-pass MVP for Q1 (commitment timing in human
fetal lung). Per notes/q1_mvp_design.md.

Pipeline:
  1. Fetch He 2022 atlas (dataset 3dc61ca1) from CELLxGENE Census
  2. Filter to 9 epithelial cell types
  3. Normalise (counts -> CP10K -> log1p)
  4. Threshold SOX2 / SOX9 at relaxed (>0), moderate (>0.5), stringent (>1.0)
  5. Classify cells per the hybrid annotation + expression rule
  6. Aggregate per stage; write q1_per_stage_fractions.csv (wide) +
     q1_long_format.csv (figure-ready, threshold-labelled)
  7. Per stage: sc.tl.rank_genes_groups between annotated_airway and
     annotated_alveolar; write q1_markers_by_stage.csv

Stops at MVP. No plotting, no Q2, no pseudotime.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "metadata"

# --- Substrate constants (locked in notes/q1_mvp_design.md) -----------------

DATASET_ID = "3dc61ca1-ce40-46b6-8337-f27260fd9a03"
VALUE_FILTER = (
    f"tissue_general == 'lung' and "
    f"dataset_id == '{DATASET_ID}'"
)

CANDIDATE_POOL = "epithelial cell of lung"

ANNOTATED_AIRWAY = {
    "lung multiciliated epithelial cell",
    "respiratory tract epithelial cell",
    "basal cell",
    "club cell",
    "lung secretory cell",
    "squamous epithelial cell",
}

ANNOTATED_ALVEOLAR = {
    "pulmonary alveolar type 1 cell",
    "pulmonary alveolar type 2 cell",
}

EPITHELIAL_TYPES = ANNOTATED_AIRWAY | ANNOTATED_ALVEOLAR | {CANDIDATE_POOL}

THRESHOLDS = {
    "relaxed": 0.0,
    "moderate": 0.5,
    "stringent": 1.0,
}
DEFAULT_THRESHOLD = "moderate"

MIN_GROUP_SIZE_FOR_MARKERS = 10
TOP_N_MARKERS = 20

STAGE_ORDER = [
    "9th week post-fertilization stage",
    "11th week post-fertilization stage",
    "15th week post-fertilization stage",
    "18th week post-fertilization stage",
    "20th week post-fertilization stage",
    "22nd week post-fertilization stage",
]


# --- Helpers ----------------------------------------------------------------

def section(title: str) -> None:
    print(f"\n=== {title} ===")


def fetch_anndata() -> ad.AnnData:
    section("Fetching He 2022 atlas from CELLxGENE Census")
    print(f"  filter: {VALUE_FILTER}")
    with cellxgene_census.open_soma(census_version="stable") as census:
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=VALUE_FILTER,
            obs_column_names=[
                "dataset_id", "development_stage", "cell_type",
                "donor_id", "is_primary_data", "assay",
            ],
        )
    print(f"  fetched: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  X dtype: {type(adata.X).__name__}")
    return adata


def filter_to_epithelial(adata: ad.AnnData) -> ad.AnnData:
    section("Filtering to 9 epithelial cell types")
    keep = adata.obs["cell_type"].astype(str).isin(EPITHELIAL_TYPES)
    n_before = adata.n_obs
    sub = adata[keep].copy()
    print(f"  {sub.n_obs:,} of {n_before:,} cells retained")
    print(f"  cell-type counts:")
    print(sub.obs["cell_type"].astype(str).value_counts().to_string())
    return sub


def normalise(adata: ad.AnnData) -> None:
    section("Normalising (CP10K + log1p)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("  done; X is now log-normalised")


def get_gene_expression(adata: ad.AnnData, gene_symbol: str) -> np.ndarray:
    """Return per-cell expression for a gene by symbol (var.feature_name)."""
    if "feature_name" in adata.var.columns:
        matches = adata.var.index[adata.var["feature_name"] == gene_symbol]
    else:
        matches = adata.var.index[adata.var.index == gene_symbol]
    if len(matches) == 0:
        raise KeyError(f"gene {gene_symbol!r} not found in adata.var")
    if len(matches) > 1:
        print(f"  warning: {len(matches)} matches for {gene_symbol!r}; using first")
    idx = adata.var.index.get_loc(matches[0])
    col = adata.X[:, idx]
    if hasattr(col, "toarray"):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col


def classify_cells(adata: ad.AnnData,
                   sox2: np.ndarray, sox9: np.ndarray) -> pd.DataFrame:
    """Assign each cell a class per the design. Returns a per-cell DataFrame
    with one column per (class, threshold) combination, plus the
    annotation-only columns."""
    cell_type = adata.obs["cell_type"].astype(str).values
    rows = pd.DataFrame(index=adata.obs.index)
    rows["cell_type"] = cell_type
    rows["stage"] = adata.obs["development_stage"].astype(str).values
    rows["donor_id"] = adata.obs["donor_id"].astype(str).values
    rows["sox2_lognorm"] = sox2
    rows["sox9_lognorm"] = sox9

    # Threshold-independent annotated classes
    rows["annotated_airway"] = np.isin(cell_type, list(ANNOTATED_AIRWAY))
    rows["annotated_alveolar"] = np.isin(cell_type, list(ANNOTATED_ALVEOLAR))
    rows["candidate_pool"] = cell_type == CANDIDATE_POOL

    # Threshold-dependent expression classes (within candidate pool only)
    for label, t in THRESHOLDS.items():
        sox2_pos = sox2 > t
        sox9_pos = sox9 > t
        in_pool = rows["candidate_pool"].values
        rows[f"candidate_bipotent_{label}"] = in_pool & sox2_pos & sox9_pos
        rows[f"transitioning_airway_{label}"] = in_pool & sox2_pos & (~sox9_pos)
        rows[f"transitioning_alveolar_{label}"] = in_pool & (~sox2_pos) & sox9_pos
        rows[f"undefined_{label}"] = in_pool & (~sox2_pos) & (~sox9_pos)

    return rows


def expression_distribution_summary(sox2: np.ndarray, sox9: np.ndarray,
                                    candidate_mask: np.ndarray) -> None:
    section("SOX2 / SOX9 distribution (within candidate pool — sanity check)")
    s2 = sox2[candidate_mask]
    s9 = sox9[candidate_mask]
    for name, vals in [("SOX2", s2), ("SOX9", s9)]:
        nz = vals[vals > 0]
        q = np.quantile(vals, [0.0, 0.25, 0.5, 0.75, 0.95, 1.0])
        print(f"  {name}: n={len(vals)}, frac_nonzero={(vals > 0).mean():.3f}, "
              f"quantiles[0,25,50,75,95,100]={q.round(3).tolist()}")
        for thr_name, thr in THRESHOLDS.items():
            print(f"    fraction > {thr} ({thr_name}): {(vals > thr).mean():.3f}")


def per_stage_table(rows: pd.DataFrame) -> pd.DataFrame:
    """Wide format: one row per stage."""
    out_records = []
    for stage in STAGE_ORDER:
        sub = rows[rows["stage"] == stage]
        n_total = len(sub)
        if n_total == 0:
            continue
        record = {
            "stage": stage,
            "n_total_epithelial": n_total,
            "n_annotated_airway": int(sub["annotated_airway"].sum()),
            "n_annotated_alveolar": int(sub["annotated_alveolar"].sum()),
            "n_candidate_pool_total": int(sub["candidate_pool"].sum()),
        }
        for label in THRESHOLDS:
            record[f"n_candidate_bipotent_{label}"] = int(sub[f"candidate_bipotent_{label}"].sum())
            record[f"n_transitioning_airway_{label}"] = int(sub[f"transitioning_airway_{label}"].sum())
            record[f"n_transitioning_alveolar_{label}"] = int(sub[f"transitioning_alveolar_{label}"].sum())
            record[f"n_undefined_{label}"] = int(sub[f"undefined_{label}"].sum())

        # Fractions (denominator = n_total_epithelial)
        for k in list(record.keys()):
            if k.startswith("n_") and k != "n_total_epithelial":
                record[f"frac_{k[2:]}"] = record[k] / n_total
        out_records.append(record)
    return pd.DataFrame(out_records)


def long_format_table(rows: pd.DataFrame) -> pd.DataFrame:
    """Long format: stage × class × threshold."""
    out_records = []
    for stage in STAGE_ORDER:
        sub = rows[rows["stage"] == stage]
        n_total = len(sub)
        if n_total == 0:
            continue
        # Threshold-independent classes
        for class_name, mask_col in [
            ("annotated_airway", "annotated_airway"),
            ("annotated_alveolar", "annotated_alveolar"),
            ("candidate_pool_total", "candidate_pool"),
        ]:
            n = int(sub[mask_col].sum())
            out_records.append({
                "stage": stage, "class": class_name, "threshold": "n/a",
                "n_cells": n, "fraction": n / n_total,
            })
        # Threshold-dependent classes
        for label in THRESHOLDS:
            for class_root in (
                "candidate_bipotent",
                "transitioning_airway",
                "transitioning_alveolar",
                "undefined",
            ):
                col = f"{class_root}_{label}"
                n = int(sub[col].sum())
                out_records.append({
                    "stage": stage, "class": class_root, "threshold": label,
                    "n_cells": n, "fraction": n / n_total,
                })
    return pd.DataFrame(out_records)


def markers_per_stage(adata: ad.AnnData, rows: pd.DataFrame) -> pd.DataFrame:
    section("Per-stage marker analysis (annotated_airway vs annotated_alveolar)")
    out_records = []
    adata = adata.copy()
    # Add the rows back to adata.obs for grouping
    adata.obs["stage_str"] = rows["stage"].values
    adata.obs["airway_alveolar_label"] = "other"
    adata.obs.loc[rows["annotated_airway"].values, "airway_alveolar_label"] = "airway"
    adata.obs.loc[rows["annotated_alveolar"].values, "airway_alveolar_label"] = "alveolar"

    for stage in STAGE_ORDER:
        stage_mask = (adata.obs["stage_str"] == stage).values
        adata_stage = adata[stage_mask].copy()
        n_air = (adata_stage.obs["airway_alveolar_label"] == "airway").sum()
        n_alv = (adata_stage.obs["airway_alveolar_label"] == "alveolar").sum()
        if n_air < MIN_GROUP_SIZE_FOR_MARKERS or n_alv < MIN_GROUP_SIZE_FOR_MARKERS:
            print(f"  {stage}: SKIP (airway={n_air}, alveolar={n_alv}; "
                  f"min={MIN_GROUP_SIZE_FOR_MARKERS})")
            continue
        # Run rank_genes_groups
        try:
            sc.tl.rank_genes_groups(
                adata_stage,
                groupby="airway_alveolar_label",
                groups=["airway", "alveolar"],
                reference="rest",
                method="wilcoxon",
                use_raw=False,
            )
        except Exception as e:
            print(f"  {stage}: rank_genes_groups failed: {e}")
            continue

        for direction in ("airway", "alveolar"):
            df = sc.get.rank_genes_groups_df(adata_stage, group=direction)
            df = df.head(TOP_N_MARKERS)
            label = "airway_top" if direction == "airway" else "alveolar_top"
            for rank_i, (_, r) in enumerate(df.iterrows(), start=1):
                gene = r["names"]
                # Map back to gene symbol if var has feature_name
                if "feature_name" in adata_stage.var.columns and gene in adata_stage.var.index:
                    sym = adata_stage.var.loc[gene, "feature_name"]
                else:
                    sym = gene
                out_records.append({
                    "stage": stage,
                    "direction": label,
                    "rank": rank_i,
                    "gene": sym,
                    "log_fold_change": float(r.get("logfoldchanges", float("nan"))),
                    "adjusted_p_value": float(r.get("pvals_adj", float("nan"))),
                })
        print(f"  {stage}: airway={n_air}, alveolar={n_alv}, "
              f"top markers extracted")
    return pd.DataFrame(out_records)


# --- Main -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tipcommit-q1-mvp")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    adata = fetch_anndata()
    adata_epi = filter_to_epithelial(adata)
    normalise(adata_epi)

    section("Per-cell SOX2 / SOX9 expression")
    sox2 = get_gene_expression(adata_epi, "SOX2")
    sox9 = get_gene_expression(adata_epi, "SOX9")
    print(f"  SOX2: nonzero in {(sox2 > 0).sum()} of {len(sox2)} cells")
    print(f"  SOX9: nonzero in {(sox9 > 0).sum()} of {len(sox9)} cells")

    rows = classify_cells(adata_epi, sox2, sox9)

    expression_distribution_summary(sox2, sox9, rows["candidate_pool"].values)

    # Per-stage tables
    section("Aggregating per stage")
    wide = per_stage_table(rows)
    long = long_format_table(rows)
    print(f"  wide format: {len(wide)} stages")
    print(f"  long format: {len(long)} (stage, class, threshold) rows")

    # Markers
    markers = markers_per_stage(adata_epi, rows)
    print(f"  markers: {len(markers)} rows total")

    # Write outputs
    wide_path = args.out_dir / "q1_per_stage_fractions.csv"
    long_path = args.out_dir / "q1_long_format.csv"
    markers_path = args.out_dir / "q1_markers_by_stage.csv"
    wide.to_csv(wide_path, index=False)
    long.to_csv(long_path, index=False)
    markers.to_csv(markers_path, index=False)

    section("Outputs written")
    for p in (wide_path, long_path, markers_path):
        try:
            print(f"  {p.relative_to(REPO_ROOT)}  ({p.stat().st_size:,} bytes)")
        except ValueError:
            print(f"  {p}  ({p.stat().st_size:,} bytes)")

    section(f"Headline: candidate_bipotent fractions at moderate threshold (default = '{DEFAULT_THRESHOLD}', > {THRESHOLDS[DEFAULT_THRESHOLD]})")
    headline = wide[["stage", "n_total_epithelial",
                     f"n_candidate_bipotent_{DEFAULT_THRESHOLD}",
                     f"frac_candidate_bipotent_{DEFAULT_THRESHOLD}"]]
    print(headline.to_string(index=False))

    section("Sensitivity: candidate_bipotent fractions at all thresholds")
    sens = wide[["stage"] + [f"frac_candidate_bipotent_{t}" for t in THRESHOLDS]]
    print(sens.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
