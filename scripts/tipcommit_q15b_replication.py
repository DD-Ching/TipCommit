#!/usr/bin/env python3
"""tipcommit q1.5b - cross-dataset replication of the Q1 v0 result on
Cao et al. 2020 (1M-cell subset). Per notes/q15b_replication_design.md.

Near-clone of scripts/tipcommit_q1_mvp.py. Differences:
  - Cao 2020 dataset_id (`fa27492b...`) instead of He 2022
  - Cao 2020 lung has no AT1/AT2 annotations; ANNOTATED_ALVEOLAR is empty
  - candidate-pool annotation is `epithelial cell of lower respiratory tract`
  - Stage range is 12-17 wpc (5 stages); 15 wpc is the only direct overlap
    with He 2022
  - Marker DE switches from annotated_airway-vs-annotated_alveolar to
    transitioning_airway-vs-transitioning_alveolar (within candidate
    pool, at moderate threshold). SOX2/SOX9 separate the groups by
    construction; flag them in the curated marker table.
  - Wide + long outputs include both within-total and within-pool
    fractions (since annotated_alveolar = 0 always, within-total
    fractions are inflated relative to He; within-pool is the
    biologically comparable denominator).
  - Curated gene_class column attached inline (no separate refinement
    script for this pass).

No new dependencies. No plotting. No Q2.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "metadata"

# --- Substrate constants (locked in notes/q15b_replication_design.md) -------

DATASET_ID = "fa27492b-82ff-4ab7-ac61-0e2b184eee67"  # Cao 2020 1M subset
VALUE_FILTER = (
    f"tissue_general == 'lung' and "
    f"dataset_id == '{DATASET_ID}'"
)

CANDIDATE_POOL = "epithelial cell of lower respiratory tract"

# Cao 2020's lung subset only annotates ciliated + squamous as committed
# airway-side epithelial cells. No basal/club/secretory annotations exist.
ANNOTATED_AIRWAY = {
    "ciliated epithelial cell",
    "squamous epithelial cell",
}

# Cao 2020 has no AT1 / AT2 annotations in lung. Empty by design.
ANNOTATED_ALVEOLAR: set[str] = set()

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
    "12th week post-fertilization stage",
    "13th week post-fertilization stage",
    "15th week post-fertilization stage",
    "16th week post-fertilization stage",
    "17th week post-fertilization stage",
]

# --- Gene-class curation (mirrors scripts/tipcommit_q1_refinement.py) -------

RIB_RE = re.compile(r"^(RPL|RPS|MRPL|MRPS)\d")
MITO_RE = re.compile(r"^MT-")
HOUSEKEEPING = {
    "MALAT1", "NEAT1", "ZFAS1", "XIST",
    "NPM1", "HNRNPA1", "NACA", "RACK1",
    "GNAS", "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSBP1",
    "MT2A", "MT1X", "MT1E", "MT1G", "MT1F",
    "TMSB4X", "TMSB10",
    "ACTB", "GAPDH", "B2M", "PPIA", "HPRT1",
}


def classify_gene(g: object) -> str:
    if not isinstance(g, str):
        return "unannotated"
    if g.startswith("ENSG"):
        return "unannotated"
    if RIB_RE.match(g):
        return "ribosomal"
    if MITO_RE.match(g):
        return "mitochondrial"
    if g in HOUSEKEEPING:
        return "housekeeping"
    return "biological"


# --- Helpers (copied verbatim from MVP where possible) ----------------------

def section(title: str) -> None:
    print(f"\n=== {title} ===")


def fetch_anndata() -> ad.AnnData:
    section("Fetching Cao 2020 (1M subset) lung cells from CELLxGENE Census")
    print(f"  filter: {VALUE_FILTER}")
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
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
    section("Filtering to Cao 2020 epithelial cell types")
    print(f"  keeping: {sorted(EPITHELIAL_TYPES)}")
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
    cell_type = adata.obs["cell_type"].astype(str).values
    rows = pd.DataFrame(index=adata.obs.index)
    rows["cell_type"] = cell_type
    rows["stage"] = adata.obs["development_stage"].astype(str).values
    rows["donor_id"] = adata.obs["donor_id"].astype(str).values
    rows["sox2_lognorm"] = sox2
    rows["sox9_lognorm"] = sox9

    rows["annotated_airway"] = np.isin(cell_type, list(ANNOTATED_AIRWAY))
    rows["annotated_alveolar"] = np.isin(cell_type, list(ANNOTATED_ALVEOLAR))
    rows["candidate_pool"] = cell_type == CANDIDATE_POOL

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
        q = np.quantile(vals, [0.0, 0.25, 0.5, 0.75, 0.95, 1.0])
        print(f"  {name}: n={len(vals)}, frac_nonzero={(vals > 0).mean():.3f}, "
              f"quantiles[0,25,50,75,95,100]={q.round(3).tolist()}")
        for thr_name, thr in THRESHOLDS.items():
            print(f"    fraction > {thr} ({thr_name}): {(vals > thr).mean():.3f}")


def per_stage_table(rows: pd.DataFrame) -> pd.DataFrame:
    """Wide format with both within-total and within-pool fractions for the
    four candidate-pool sub-states (since Cao's empty annotated_alveolar
    inflates within-total fractions relative to He 2022)."""
    out_records = []
    for stage in STAGE_ORDER:
        sub = rows[rows["stage"] == stage]
        n_total = len(sub)
        if n_total == 0:
            continue
        n_pool = int(sub["candidate_pool"].sum())
        record = {
            "stage": stage,
            "n_total_epithelial": n_total,
            "n_annotated_airway": int(sub["annotated_airway"].sum()),
            "n_annotated_alveolar": int(sub["annotated_alveolar"].sum()),
            "n_candidate_pool_total": n_pool,
        }
        for label in THRESHOLDS:
            for sub_state in (
                "candidate_bipotent",
                "transitioning_airway",
                "transitioning_alveolar",
                "undefined",
            ):
                col = f"{sub_state}_{label}"
                record[f"n_{col}"] = int(sub[col].sum())

        # Within-total-epithelial fractions
        for k in list(record.keys()):
            if k.startswith("n_") and k != "n_total_epithelial":
                record[f"frac_total_{k[2:]}"] = record[k] / n_total

        # Within-candidate-pool fractions for the four sub-states
        if n_pool > 0:
            for label in THRESHOLDS:
                for sub_state in (
                    "candidate_bipotent",
                    "transitioning_airway",
                    "transitioning_alveolar",
                    "undefined",
                ):
                    record[f"frac_pool_{sub_state}_{label}"] = (
                        record[f"n_{sub_state}_{label}"] / n_pool
                    )
        out_records.append(record)
    return pd.DataFrame(out_records)


def long_format_table(rows: pd.DataFrame) -> pd.DataFrame:
    """Long format with explicit denominator + threshold columns."""
    out_records = []
    for stage in STAGE_ORDER:
        sub = rows[rows["stage"] == stage]
        n_total = len(sub)
        if n_total == 0:
            continue
        n_pool = int(sub["candidate_pool"].sum())

        # Threshold-independent classes (within-total only)
        for class_name, mask_col in [
            ("annotated_airway", "annotated_airway"),
            ("annotated_alveolar", "annotated_alveolar"),
            ("candidate_pool_total", "candidate_pool"),
        ]:
            n = int(sub[mask_col].sum())
            out_records.append({
                "stage": stage, "class": class_name, "threshold": "n/a",
                "denominator": "total_epithelial",
                "n_cells": n, "fraction": n / n_total,
            })

        # Threshold-dependent classes: both denominators
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
                    "denominator": "total_epithelial",
                    "n_cells": n, "fraction": n / n_total,
                })
                out_records.append({
                    "stage": stage, "class": class_root, "threshold": label,
                    "denominator": "candidate_pool",
                    "n_cells": n,
                    "fraction": (n / n_pool) if n_pool > 0 else float("nan"),
                })
    return pd.DataFrame(out_records)


def markers_per_stage(adata: ad.AnnData, rows: pd.DataFrame) -> pd.DataFrame:
    """Replication marker DE: transitioning_airway vs transitioning_alveolar
    within the candidate pool, defined at the moderate threshold.

    SOX2 and SOX9 themselves separate the groups by construction; the
    curated table flags them as such (see classify_gene + the marker
    note in notes/q15b_replication_design.md)."""
    section(
        "Per-stage marker analysis: transitioning_airway vs "
        "transitioning_alveolar (within candidate pool, moderate threshold)"
    )
    out_records = []
    adata = adata.copy()
    adata.obs["stage_str"] = rows["stage"].values
    label_arr = np.full(adata.n_obs, "other", dtype=object)
    label_arr[rows[f"transitioning_airway_{DEFAULT_THRESHOLD}"].values] = "airway"
    label_arr[rows[f"transitioning_alveolar_{DEFAULT_THRESHOLD}"].values] = "alveolar"
    adata.obs["airway_alveolar_label"] = label_arr

    for stage in STAGE_ORDER:
        stage_mask = (adata.obs["stage_str"] == stage).values
        adata_stage = adata[stage_mask].copy()
        n_air = (adata_stage.obs["airway_alveolar_label"] == "airway").sum()
        n_alv = (adata_stage.obs["airway_alveolar_label"] == "alveolar").sum()
        if n_air < MIN_GROUP_SIZE_FOR_MARKERS or n_alv < MIN_GROUP_SIZE_FOR_MARKERS:
            print(f"  {stage}: SKIP (airway={n_air}, alveolar={n_alv}; "
                  f"min={MIN_GROUP_SIZE_FOR_MARKERS})")
            continue
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

    df_out = pd.DataFrame(out_records)
    if len(df_out) > 0:
        df_out["gene_class"] = df_out["gene"].map(classify_gene)
        df_out["interpretable"] = df_out["gene_class"] == "biological"
        # Recurrence count (biological only) per direction
        rec = (
            df_out[df_out["interpretable"]]
            .groupby(["direction", "gene"])["stage"]
            .nunique()
            .rename("recurrence_n_stages")
            .reset_index()
        )
        df_out = df_out.merge(rec, on=["direction", "gene"], how="left")
        df_out["recurrence_n_stages"] = df_out["recurrence_n_stages"].fillna(0).astype(int)
    return df_out


# --- Main -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tipcommit-q15b-replication")
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

    section("Aggregating per stage")
    wide = per_stage_table(rows)
    long = long_format_table(rows)
    print(f"  wide format: {len(wide)} stages")
    print(f"  long format: {len(long)} (stage, class, threshold, denom) rows")

    markers = markers_per_stage(adata_epi, rows)
    print(f"  markers: {len(markers)} rows total")

    wide_path = args.out_dir / "q15b_replication_per_stage_fractions.csv"
    long_path = args.out_dir / "q15b_replication_long_format.csv"
    markers_path = args.out_dir / "q15b_replication_markers_by_stage.csv"
    wide.to_csv(wide_path, index=False)
    long.to_csv(long_path, index=False)
    markers.to_csv(markers_path, index=False)

    section("Outputs written")
    for p in (wide_path, long_path, markers_path):
        try:
            print(f"  {p.relative_to(REPO_ROOT)}  ({p.stat().st_size:,} bytes)")
        except ValueError:
            print(f"  {p}  ({p.stat().st_size:,} bytes)")

    section(f"Headline: candidate_bipotent at moderate threshold (>{THRESHOLDS[DEFAULT_THRESHOLD]})")
    cols_total = [
        "stage", "n_total_epithelial", "n_candidate_pool_total",
        f"n_candidate_bipotent_{DEFAULT_THRESHOLD}",
        f"frac_total_candidate_bipotent_{DEFAULT_THRESHOLD}",
        f"frac_pool_candidate_bipotent_{DEFAULT_THRESHOLD}",
    ]
    print(wide[cols_total].to_string(index=False))

    section("Sensitivity: candidate_bipotent within-pool fraction at all thresholds")
    cols_sens = ["stage"] + [f"frac_pool_candidate_bipotent_{t}" for t in THRESHOLDS]
    print(wide[cols_sens].to_string(index=False))

    if len(markers) > 0:
        section("Curated marker recurrence (biological only, >=2 stages)")
        bio = markers[markers["interpretable"]]
        for direction in ("airway_top", "alveolar_top"):
            print(f"\n  --- {direction} ---")
            sub = bio[bio["direction"] == direction]
            grouped = (
                sub.groupby("gene", observed=True)
                .agg(
                    n_stages=("stage", "nunique"),
                    min_rank=("rank", "min"),
                )
                .query("n_stages >= 2")
                .sort_values(["n_stages", "min_rank"], ascending=[False, True])
            )
            print(grouped.head(20).to_string())

    return 0


if __name__ == "__main__":
    sys.exit(main())
