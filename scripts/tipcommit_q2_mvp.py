#!/usr/bin/env python3
"""tipcommit q2-mvp - axis continuity (gradient vs switch vs mixed).

Per notes/q2_design.md, with the implementation refinement that the
commitment diagnostics are evaluated in TWO views per substrate:
  1. all_epithelial - the full epithelial cell set
  2. exiting_tip   - cells with distal_tip_score below the substrate-wide
                     median (so we don't overinterpret commitment in cells
                     still strongly tip-like)

Pipeline:
  1. Fetch He 2022 epithelial slice (primary)
  2. Fetch Cao 2020 epithelial slice (secondary, 15-17 wpc only)
  3. Normalise each independently
  4. Score each cell on the three locked signatures via
     sc.tl.score_genes
  5. Per (substrate, view, stage), compute:
       - Sarle bimodality coefficient on commitment = airway - alveolar
       - gap-zone fraction
       - co-commitment fraction
  6. 5x5 binned joint density of (airway_score, alveolar_score) per
     (substrate, view, stage)
  7. Write 4 CSVs

No plotting. No clustering. No pseudotime. No new dependencies.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import skew, kurtosis

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "metadata"

# --- Locked signatures (notes/q2_design.md section C) -----------------------

DISTAL_TIP_GENES = ["SOX9", "ID2", "ETV5", "TESC"]
AIRWAY_GENES = ["SOX2", "TP63", "SCGB3A2", "FOXJ1"]
ALVEOLAR_GENES = ["SFTPC", "SFTPB", "AGER", "SLC34A2"]
ALL_SIG_GENES = DISTAL_TIP_GENES + AIRWAY_GENES + ALVEOLAR_GENES

# --- Substrate constants ----------------------------------------------------

HE_DATASET_ID = "3dc61ca1-ce40-46b6-8337-f27260fd9a03"
HE_VALUE_FILTER = (
    f"tissue_general == 'lung' and dataset_id == '{HE_DATASET_ID}'"
)
HE_EPITHELIAL_TYPES = {
    "epithelial cell of lung",
    "lung multiciliated epithelial cell",
    "respiratory tract epithelial cell",
    "basal cell",
    "club cell",
    "lung secretory cell",
    "squamous epithelial cell",
    "pulmonary alveolar type 1 cell",
    "pulmonary alveolar type 2 cell",
}
HE_STAGE_ORDER = [
    "9th week post-fertilization stage",
    "11th week post-fertilization stage",
    "15th week post-fertilization stage",
    "18th week post-fertilization stage",
    "20th week post-fertilization stage",
    "22nd week post-fertilization stage",
]

CAO_DATASET_ID = "fa27492b-82ff-4ab7-ac61-0e2b184eee67"
CAO_VALUE_FILTER = (
    f"tissue_general == 'lung' and dataset_id == '{CAO_DATASET_ID}'"
)
CAO_EPITHELIAL_TYPES = {
    "epithelial cell of lower respiratory tract",
    "ciliated epithelial cell",
    "squamous epithelial cell",
}
CAO_STAGE_ORDER_KEEP = [
    "15th week post-fertilization stage",
    "16th week post-fertilization stage",
    "17th week post-fertilization stage",
]

# --- Diagnostic constants ---------------------------------------------------

GAP_ZONE_COMMIT_FRAC = 0.5    # |commitment| <= GAP_ZONE_COMMIT_FRAC * max
GAP_ZONE_DISTAL_PCT = 0.75    # distal_tip_score below this percentile
N_DENSITY_BINS = 5            # 5x5 joint density grid
EXITING_TIP_PCT = 0.5         # below median = exiting-tip view


# --- Helpers (mirrors q1_mvp + q15b helpers) --------------------------------

def section(title: str) -> None:
    print(f"\n=== {title} ===")


def fetch_anndata(value_filter: str, label: str) -> ad.AnnData:
    section(f"Fetching {label} from CELLxGENE Census")
    print(f"  filter: {value_filter}")
    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=value_filter,
            obs_column_names=[
                "dataset_id", "development_stage", "cell_type",
                "donor_id", "is_primary_data", "assay",
            ],
        )
    print(f"  fetched: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def filter_to_epithelial(adata: ad.AnnData, types: set[str], label: str) -> ad.AnnData:
    section(f"{label}: filtering to epithelial cell types")
    keep = adata.obs["cell_type"].astype(str).isin(types)
    sub = adata[keep].copy()
    print(f"  {sub.n_obs:,} of {adata.n_obs:,} cells retained")
    print(sub.obs["cell_type"].astype(str).value_counts().to_string())
    return sub


def normalise(adata: ad.AnnData) -> None:
    section("Normalising (CP10K + log1p)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def map_symbols_to_var_index(adata: ad.AnnData,
                             symbols: list[str]) -> dict[str, str]:
    """Return {gene_symbol: adata.var index} for each symbol; raise if any
    symbol is missing (per design's 'fail loudly' rule)."""
    if "feature_name" not in adata.var.columns:
        raise RuntimeError("adata.var has no feature_name column")
    out: dict[str, str] = {}
    missing: list[str] = []
    for sym in symbols:
        matches = adata.var.index[adata.var["feature_name"] == sym]
        if len(matches) == 0:
            missing.append(sym)
        else:
            out[sym] = matches[0]
    if missing:
        raise KeyError(
            f"signature genes not found in adata.var: {missing}"
        )
    return out


def score_signatures(adata: ad.AnnData, label: str) -> None:
    section(f"{label}: scoring three locked signatures")
    sym_to_var = map_symbols_to_var_index(adata, ALL_SIG_GENES)
    for score_name, gene_list in [
        ("distal_tip_score", DISTAL_TIP_GENES),
        ("airway_score", AIRWAY_GENES),
        ("alveolar_score", ALVEOLAR_GENES),
    ]:
        var_ids = [sym_to_var[g] for g in gene_list]
        sc.tl.score_genes(
            adata,
            gene_list=var_ids,
            score_name=score_name,
            n_bins=25,
            ctrl_size=50,
            random_state=0,
        )
        print(f"  {score_name}: genes={gene_list}, "
              f"range=[{adata.obs[score_name].min():.3f}, "
              f"{adata.obs[score_name].max():.3f}]")
    adata.obs["commitment_score"] = (
        adata.obs["airway_score"] - adata.obs["alveolar_score"]
    )


# --- Diagnostic functions ---------------------------------------------------

def bimodality_coefficient(x: np.ndarray) -> float:
    """Sarle's bimodality coefficient. BC > 5/9 (0.555...) is conventional
    evidence of bimodality."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    n = len(x)
    if n < 4:
        return float("nan")
    g = float(skew(x))
    k = float(kurtosis(x, fisher=True))
    denom = k + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    if denom == 0:
        return float("nan")
    return (g ** 2 + 1.0) / denom


def diagnostics_for(commit: np.ndarray,
                    distal: np.ndarray,
                    airway: np.ndarray,
                    alveolar: np.ndarray,
                    commit_max_abs: float,
                    distal_q75: float,
                    airway_med: float,
                    alveolar_med: float) -> dict:
    """Three diagnostics + counts on the cells passed in.

    Reference quantities (commit_max_abs, distal_q75, airway_med,
    alveolar_med) are computed substrate-wide so per-stage diagnostics
    are comparable across stages."""
    n = len(commit)
    if n < 4:
        return dict(n_cells=n, bimodality_coefficient=float("nan"),
                    gap_zone_fraction=float("nan"),
                    co_commitment_fraction=float("nan"))
    bc = bimodality_coefficient(commit)
    gap = ((np.abs(commit) <= GAP_ZONE_COMMIT_FRAC * commit_max_abs)
           & (distal < distal_q75))
    co = (airway > airway_med) & (alveolar > alveolar_med)
    return dict(
        n_cells=n,
        bimodality_coefficient=bc,
        gap_zone_fraction=float(gap.mean()),
        co_commitment_fraction=float(co.mean()),
    )


def per_stage_diagnostics(obs: pd.DataFrame,
                          stages: list[str],
                          substrate_label: str,
                          exiting_tip_threshold: float,
                          commit_max_abs: float,
                          distal_q75: float,
                          airway_med: float,
                          alveolar_med: float) -> pd.DataFrame:
    rows = []
    obs = obs.copy()
    obs["in_exiting_tip"] = obs["distal_tip_score"] < exiting_tip_threshold
    for stage in stages:
        stage_mask = (obs["development_stage"].astype(str) == stage).values
        if not stage_mask.any():
            continue
        for view, mask in [
            ("all_epithelial", stage_mask),
            ("exiting_tip", stage_mask & obs["in_exiting_tip"].values),
        ]:
            sub = obs[mask]
            d = diagnostics_for(
                commit=sub["commitment_score"].values,
                distal=sub["distal_tip_score"].values,
                airway=sub["airway_score"].values,
                alveolar=sub["alveolar_score"].values,
                commit_max_abs=commit_max_abs,
                distal_q75=distal_q75,
                airway_med=airway_med,
                alveolar_med=alveolar_med,
            )
            rows.append({
                "substrate": substrate_label,
                "view": view,
                "stage": stage,
                **d,
                "median_distal_tip_score": float(sub["distal_tip_score"].median()) if len(sub) else float("nan"),
                "median_airway_score": float(sub["airway_score"].median()) if len(sub) else float("nan"),
                "median_alveolar_score": float(sub["alveolar_score"].median()) if len(sub) else float("nan"),
            })
    return pd.DataFrame(rows)


def joint_density_grid(obs: pd.DataFrame,
                       stages: list[str],
                       substrate_label: str,
                       airway_edges: np.ndarray,
                       alveolar_edges: np.ndarray,
                       exiting_tip_threshold: float) -> pd.DataFrame:
    """5x5 binned joint density of (airway_score, alveolar_score) per
    (substrate, view, stage). Bin edges are substrate-wide so bins are
    comparable across stages."""
    rows = []
    obs = obs.copy()
    obs["in_exiting_tip"] = obs["distal_tip_score"] < exiting_tip_threshold
    for stage in stages:
        stage_mask = (obs["development_stage"].astype(str) == stage).values
        if not stage_mask.any():
            continue
        for view, mask in [
            ("all_epithelial", stage_mask),
            ("exiting_tip", stage_mask & obs["in_exiting_tip"].values),
        ]:
            sub = obs[mask]
            n = len(sub)
            if n == 0:
                continue
            ai = np.clip(np.digitize(sub["airway_score"], airway_edges[1:-1]),
                         0, N_DENSITY_BINS - 1)
            av = np.clip(np.digitize(sub["alveolar_score"], alveolar_edges[1:-1]),
                         0, N_DENSITY_BINS - 1)
            cnt = pd.crosstab(pd.Series(ai, name="airway_bin"),
                              pd.Series(av, name="alveolar_bin")).reindex(
                index=range(N_DENSITY_BINS), columns=range(N_DENSITY_BINS),
                fill_value=0,
            )
            for ai_b in range(N_DENSITY_BINS):
                for av_b in range(N_DENSITY_BINS):
                    rows.append({
                        "substrate": substrate_label,
                        "view": view,
                        "stage": stage,
                        "airway_bin": ai_b,
                        "alveolar_bin": av_b,
                        "n_cells": int(cnt.loc[ai_b, av_b]),
                        "fraction_of_view": float(cnt.loc[ai_b, av_b]) / n,
                    })
    return pd.DataFrame(rows)


@dataclass
class SubstrateOutputs:
    label: str
    obs: pd.DataFrame
    diagnostics: pd.DataFrame
    density: pd.DataFrame
    exiting_tip_threshold: float


def process_substrate(value_filter: str,
                      types: set[str],
                      stages: list[str],
                      label: str) -> SubstrateOutputs:
    adata = fetch_anndata(value_filter, label)
    adata_epi = filter_to_epithelial(adata, types, label)
    normalise(adata_epi)
    score_signatures(adata_epi, label)

    obs = adata_epi.obs[[
        "development_stage", "cell_type", "donor_id", "assay",
        "distal_tip_score", "airway_score", "alveolar_score",
        "commitment_score",
    ]].copy()
    obs["substrate"] = label
    obs["development_stage"] = obs["development_stage"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)

    # Substrate-wide reference quantities
    distal_med = float(obs["distal_tip_score"].median())
    distal_q75 = float(obs["distal_tip_score"].quantile(0.75))
    airway_med = float(obs["airway_score"].median())
    alveolar_med = float(obs["alveolar_score"].median())
    commit_max_abs = float(np.abs(obs["commitment_score"]).max())
    section(f"{label}: substrate-wide reference quantities")
    print(f"  distal_tip_score median: {distal_med:.4f}  (exiting_tip cutoff)")
    print(f"  distal_tip_score q75:    {distal_q75:.4f}  (gap-zone cutoff)")
    print(f"  airway_score median:     {airway_med:.4f}")
    print(f"  alveolar_score median:   {alveolar_med:.4f}")
    print(f"  |commitment_score| max:  {commit_max_abs:.4f}")

    diag = per_stage_diagnostics(
        obs=obs,
        stages=stages,
        substrate_label=label,
        exiting_tip_threshold=distal_med,
        commit_max_abs=commit_max_abs,
        distal_q75=distal_q75,
        airway_med=airway_med,
        alveolar_med=alveolar_med,
    )

    # Substrate-wide bin edges for the joint-density grid (quintiles).
    airway_edges = np.quantile(
        obs["airway_score"], np.linspace(0, 1, N_DENSITY_BINS + 1)
    )
    alveolar_edges = np.quantile(
        obs["alveolar_score"], np.linspace(0, 1, N_DENSITY_BINS + 1)
    )
    density = joint_density_grid(
        obs=obs,
        stages=stages,
        substrate_label=label,
        airway_edges=airway_edges,
        alveolar_edges=alveolar_edges,
        exiting_tip_threshold=distal_med,
    )

    obs_out = obs.copy()
    obs_out["in_exiting_tip"] = obs_out["distal_tip_score"] < distal_med
    return SubstrateOutputs(
        label=label, obs=obs_out, diagnostics=diag, density=density,
        exiting_tip_threshold=distal_med,
    )


# --- Verdict scoring (pre-registered, per notes/q2_design.md section E) -----

BC_THRESHOLD = 5.0 / 9.0
GAP_SWITCH_MAX = 0.10
CO_SWITCH_MAX = 0.15
GAP_GRADIENT_MIN = 0.25
CO_GRADIENT_MIN = 0.30


def score_verdict(diag: pd.DataFrame, view: str) -> dict:
    """Apply the design's pre-registered verdict criteria to a single
    (substrate, view) slice."""
    sub = diag[diag["view"] == view].copy()
    sub = sub.dropna(subset=["bimodality_coefficient"])
    n_stages = len(sub)
    bc_high = (sub["bimodality_coefficient"] >= BC_THRESHOLD).sum()
    bc_low = (sub["bimodality_coefficient"] < BC_THRESHOLD).sum()
    gap_low = (sub["gap_zone_fraction"] <= GAP_SWITCH_MAX).sum()
    gap_high = (sub["gap_zone_fraction"] >= GAP_GRADIENT_MIN).sum()
    co_low = (sub["co_commitment_fraction"] <= CO_SWITCH_MAX).sum()
    co_high = (sub["co_commitment_fraction"] >= CO_GRADIENT_MIN).sum()

    switch_pass = (bc_high >= 3) and (gap_low >= 3) and (co_low >= 3)
    gradient_pass = (bc_low >= 4) and (gap_high >= 3) and (co_high >= 3)

    if switch_pass and gradient_pass:
        verdict = "ambiguous (both criteria pass)"
    elif switch_pass:
        verdict = "switch-like"
    elif gradient_pass:
        verdict = "gradient-like"
    else:
        verdict = "mixed_or_inconclusive"
    return dict(
        view=view, n_stages=n_stages,
        bc_high=int(bc_high), bc_low=int(bc_low),
        gap_low=int(gap_low), gap_high=int(gap_high),
        co_low=int(co_low), co_high=int(co_high),
        switch_pass=switch_pass, gradient_pass=gradient_pass,
        verdict=verdict,
    )


# --- Main -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tipcommit-q2-mvp")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    he = process_substrate(
        HE_VALUE_FILTER, HE_EPITHELIAL_TYPES, HE_STAGE_ORDER, "he2022",
    )
    cao = process_substrate(
        CAO_VALUE_FILTER, CAO_EPITHELIAL_TYPES, CAO_STAGE_ORDER_KEEP, "cao2020",
    )

    # Per-cell scores -> long CSV
    per_cell = pd.concat([he.obs, cao.obs], axis=0)
    per_cell.index.name = "cell_id"
    per_cell = per_cell.reset_index()

    # Diagnostics: He goes into the primary CSV; Cao into its own.
    diag_he = he.diagnostics
    diag_cao = cao.diagnostics
    density_all = pd.concat([he.density, cao.density], axis=0).reset_index(drop=True)

    per_cell_path = args.out_dir / "q2_per_cell_scores.csv"
    diag_he_path = args.out_dir / "q2_per_stage_diagnostics.csv"
    density_path = args.out_dir / "q2_joint_density.csv"
    diag_cao_path = args.out_dir / "q2_cao_secondary_diagnostics.csv"

    per_cell.to_csv(per_cell_path, index=False)
    diag_he.to_csv(diag_he_path, index=False)
    diag_cao.to_csv(diag_cao_path, index=False)
    density_all.to_csv(density_path, index=False)

    section("Outputs written")
    for p in (per_cell_path, diag_he_path, density_path, diag_cao_path):
        try:
            print(f"  {p.relative_to(REPO_ROOT)}  ({p.stat().st_size:,} bytes)")
        except ValueError:
            print(f"  {p}  ({p.stat().st_size:,} bytes)")

    # Console summaries -------------------------------------------------------
    section("He 2022 per-stage diagnostics")
    print(diag_he.round(4).to_string(index=False))

    section("Cao 2020 per-stage diagnostics (secondary; 15-17 wpc only)")
    print(diag_cao.round(4).to_string(index=False))

    section("Verdict scoring (He 2022 - primary)")
    for view in ("all_epithelial", "exiting_tip"):
        v = score_verdict(diag_he, view)
        print(f"  view={view}: verdict={v['verdict']} "
              f"(switch_pass={v['switch_pass']}, gradient_pass={v['gradient_pass']}; "
              f"bc_high={v['bc_high']}/{v['n_stages']}, "
              f"bc_low={v['bc_low']}/{v['n_stages']}; "
              f"gap_low={v['gap_low']}, gap_high={v['gap_high']}; "
              f"co_low={v['co_low']}, co_high={v['co_high']})")

    section("Verdict scoring (Cao 2020 - secondary; for direction only)")
    for view in ("all_epithelial", "exiting_tip"):
        v = score_verdict(diag_cao, view)
        print(f"  view={view}: verdict={v['verdict']} "
              f"(bc_high={v['bc_high']}/{v['n_stages']}, "
              f"bc_low={v['bc_low']}/{v['n_stages']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
