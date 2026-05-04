#!/usr/bin/env python3
"""tipcommit_q1_refinement.py - Q1.1 refinement pass.

Reads outputs of tipcommit_q1_mvp.py and produces:
- a cleaner state breakdown per stage at the moderate threshold,
  reporting fractions both within total epithelial and within
  candidate_pool only (Goal A + Goal C)
- a curated markers table flagging ribosomal / mitochondrial /
  housekeeping genes and counting per-direction recurrence (Goal B)

Pure pandas; no Census fetch; reproducible from the existing
metadata/ CSVs alone.

Outputs:
- metadata/q1_state_breakdown_by_stage.csv
- metadata/q1_markers_by_stage_curated.csv
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
META = REPO_ROOT / "metadata"

STAGE_ORDER = [
    "9th week post-fertilization stage",
    "11th week post-fertilization stage",
    "15th week post-fertilization stage",
    "18th week post-fertilization stage",
    "20th week post-fertilization stage",
    "22nd week post-fertilization stage",
]
STAGE_SHORT = {s: s.split(" ")[0] for s in STAGE_ORDER}  # "9th" / "11th" / ...

CANDIDATE_POOL_STATES = [
    "candidate_bipotent",
    "transitioning_airway",
    "transitioning_alveolar",
    "undefined",
]
ANNOTATED_STATES = ["annotated_airway", "annotated_alveolar"]
ALL_STATES = CANDIDATE_POOL_STATES + ANNOTATED_STATES

# Gene-class filters used to flag (not remove) noise in the markers table.
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


def state_breakdown() -> pd.DataFrame:
    df = pd.read_csv(META / "q1_per_stage_fractions.csv")
    rows = []
    for _, r in df.iterrows():
        n_total = r["n_total_epithelial"]
        n_pool = r["n_candidate_pool_total"]
        for state in ALL_STATES:
            if state in CANDIDATE_POOL_STATES:
                n = int(r[f"n_{state}_moderate"])
                frac_pool = n / n_pool if n_pool > 0 else float("nan")
            else:
                n = int(r[f"n_{state}"])
                frac_pool = float("nan")
            rows.append({
                "stage": r["stage"],
                "state": state,
                "n_cells": n,
                "frac_of_total_epithelial": n / n_total,
                "frac_of_candidate_pool": frac_pool,
            })
    out = pd.DataFrame(rows)
    out["stage"] = pd.Categorical(out["stage"], categories=STAGE_ORDER, ordered=True)
    out["state"] = pd.Categorical(out["state"], categories=ALL_STATES, ordered=True)
    out = out.sort_values(["stage", "state"]).reset_index(drop=True)
    return out


def curate_markers() -> pd.DataFrame:
    df = pd.read_csv(META / "q1_markers_by_stage.csv")
    df["gene_class"] = df["gene"].map(classify_gene)
    df["interpretable"] = df["gene_class"] == "biological"
    rec = (
        df[df["interpretable"]]
        .groupby(["direction", "gene"])["stage"]
        .nunique()
        .rename("recurrence_n_stages")
        .reset_index()
    )
    df = df.merge(rec, on=["direction", "gene"], how="left")
    df["recurrence_n_stages"] = df["recurrence_n_stages"].fillna(0).astype(int)
    df["stage"] = pd.Categorical(df["stage"], categories=STAGE_ORDER, ordered=True)
    df = df.sort_values(["stage", "direction", "rank"]).reset_index(drop=True)
    return df


def stages_for(stages_iter) -> str:
    s = set(stages_iter)
    return ",".join(STAGE_SHORT[x] for x in STAGE_ORDER if x in s)


def report(state_df: pd.DataFrame, marker_df: pd.DataFrame) -> None:
    print("\n=== Goal C: candidate_bipotent fraction (moderate threshold) ===")
    cb = state_df[state_df["state"] == "candidate_bipotent"][
        ["stage", "n_cells", "frac_of_total_epithelial", "frac_of_candidate_pool"]
    ].copy()
    cb["frac_of_total_epithelial"] = cb["frac_of_total_epithelial"].round(4)
    cb["frac_of_candidate_pool"] = cb["frac_of_candidate_pool"].round(4)
    print(cb.to_string(index=False))

    print("\n=== Curated airway markers (biological only, >=2 stages) ===")
    aw_df = marker_df[
        (marker_df["direction"] == "airway_top") & (marker_df["interpretable"])
    ]
    aw = (
        aw_df.groupby("gene", observed=True)
        .agg(
            n_stages=("stage", "nunique"),
            min_rank=("rank", "min"),
            stages=("stage", lambda s: stages_for(s)),
        )
        .query("n_stages >= 2")
        .sort_values(["n_stages", "min_rank"], ascending=[False, True])
    )
    print(aw.head(25).to_string())

    print("\n=== Curated alveolar markers (biological only, >=2 stages) ===")
    al_df = marker_df[
        (marker_df["direction"] == "alveolar_top") & (marker_df["interpretable"])
    ]
    al = (
        al_df.groupby("gene", observed=True)
        .agg(
            n_stages=("stage", "nunique"),
            min_rank=("rank", "min"),
            stages=("stage", lambda s: stages_for(s)),
        )
        .query("n_stages >= 2")
        .sort_values(["n_stages", "min_rank"], ascending=[False, True])
    )
    print(al.head(25).to_string())

    print("\n=== gene_class counts in markers table ===")
    print(marker_df["gene_class"].value_counts().to_string())


def main() -> int:
    state_df = state_breakdown()
    marker_df = curate_markers()
    state_df.to_csv(META / "q1_state_breakdown_by_stage.csv", index=False)
    marker_df.to_csv(META / "q1_markers_by_stage_curated.csv", index=False)
    print(f"wrote metadata/q1_state_breakdown_by_stage.csv ({len(state_df)} rows)")
    print(f"wrote metadata/q1_markers_by_stage_curated.csv ({len(marker_df)} rows)")
    report(state_df, marker_df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
