#!/usr/bin/env python3
"""tipcommit q2.1 - candidate-pool-only diagnostic pass.

Narrow follow-up to Q2 MVP. Tests the section 4a ambiguity in
notes/q2_mvp.md: is the late-stage bimodality real cell-state
hardening, or is it driven by the increasing fraction of upstream-
annotated committed cells (AT1, AT2, basal, club, multiciliated,
secretory) at later stages?

Approach: re-run the same Q2 diagnostics on the same per-cell
scores, but restricted to cell_type == 'epithelial cell of lung'
in He 2022 (the candidate pool annotation, by Q1's definition).

No Census re-fetch. No new dependencies. No new design. Reads
metadata/q2_per_cell_scores.csv produced by tipcommit_q2_mvp.py.

Substrate scope: He 2022 only. Cao 2020 was already shown
unsuitable for shape-level claims (Q1.5b).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis

REPO_ROOT = Path(__file__).resolve().parent.parent
META = REPO_ROOT / "metadata"

CANDIDATE_POOL_CT = "epithelial cell of lung"  # He 2022 annotation

HE_STAGE_ORDER = [
    "9th week post-fertilization stage",
    "11th week post-fertilization stage",
    "15th week post-fertilization stage",
    "18th week post-fertilization stage",
    "20th week post-fertilization stage",
    "22nd week post-fertilization stage",
]

# Same diagnostic constants as Q2 MVP.
GAP_ZONE_COMMIT_FRAC = 0.5
GAP_ZONE_DISTAL_PCT = 0.75
EXITING_TIP_PCT = 0.5
BC_THRESHOLD = 5.0 / 9.0


def bimodality_coefficient(x: np.ndarray) -> float:
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


def diagnostics_for(commit, distal, airway, alveolar,
                    commit_max_abs, distal_q75,
                    airway_med, alveolar_med) -> dict:
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


def main() -> int:
    in_path = META / "q2_per_cell_scores.csv"
    if not in_path.exists():
        print(f"missing input: {in_path}", file=sys.stderr)
        return 1
    print(f"reading {in_path.relative_to(REPO_ROOT)}")
    df = pd.read_csv(in_path)

    # Filter: He 2022 + candidate pool annotation only.
    pool = df[(df["substrate"] == "he2022")
              & (df["cell_type"] == CANDIDATE_POOL_CT)].copy()
    print(f"candidate pool cells: {len(pool):,} (of {len(df):,} total scored cells)")

    # Within-pool reference quantities (the candidate pool is now the universe).
    distal_med_pool = float(pool["distal_tip_score"].median())
    distal_q75_pool = float(pool["distal_tip_score"].quantile(0.75))
    airway_med_pool = float(pool["airway_score"].median())
    alveolar_med_pool = float(pool["alveolar_score"].median())
    commit_max_abs_pool = float(np.abs(pool["commitment_score"]).max())
    print(f"within-pool reference quantities:")
    print(f"  distal_tip median: {distal_med_pool:.4f}  (within-pool exiting-tip cutoff)")
    print(f"  distal_tip q75:    {distal_q75_pool:.4f}  (within-pool gap-zone cutoff)")
    print(f"  airway median:     {airway_med_pool:.4f}")
    print(f"  alveolar median:   {alveolar_med_pool:.4f}")
    print(f"  |commitment| max:  {commit_max_abs_pool:.4f}")

    pool["in_exiting_tip_pool"] = pool["distal_tip_score"] < distal_med_pool

    rows = []
    for stage in HE_STAGE_ORDER:
        stage_mask = (pool["development_stage"].astype(str) == stage).values
        if not stage_mask.any():
            continue
        for view, mask in [
            ("all_in_pool", stage_mask),
            ("exiting_tip_in_pool", stage_mask & pool["in_exiting_tip_pool"].values),
        ]:
            sub = pool[mask]
            d = diagnostics_for(
                commit=sub["commitment_score"].values,
                distal=sub["distal_tip_score"].values,
                airway=sub["airway_score"].values,
                alveolar=sub["alveolar_score"].values,
                commit_max_abs=commit_max_abs_pool,
                distal_q75=distal_q75_pool,
                airway_med=airway_med_pool,
                alveolar_med=alveolar_med_pool,
            )
            rows.append({
                "substrate": "he2022",
                "view": view,
                "stage": stage,
                **d,
                "median_distal_tip_score": float(sub["distal_tip_score"].median()) if len(sub) else float("nan"),
                "median_airway_score": float(sub["airway_score"].median()) if len(sub) else float("nan"),
                "median_alveolar_score": float(sub["alveolar_score"].median()) if len(sub) else float("nan"),
            })
    out = pd.DataFrame(rows)
    out_path = META / "q21_candidate_pool_diagnostics.csv"
    out.to_csv(out_path, index=False)
    print(f"\nwrote {out_path.relative_to(REPO_ROOT)} ({out_path.stat().st_size} bytes)")

    # Console comparison: Q2 MVP all_epithelial BC vs Q2.1 within-pool BC
    print("\n=== Per-stage diagnostics, candidate-pool only ===")
    print(out.round(4).to_string(index=False))

    # Side-by-side BC comparison with Q2 MVP
    q2_path = META / "q2_per_stage_diagnostics.csv"
    q2 = pd.read_csv(q2_path)
    print("\n=== BC side-by-side: Q2 MVP vs Q2.1 candidate-pool only (He 2022) ===")
    bc_q2_all = q2[(q2["substrate"] == "he2022") & (q2["view"] == "all_epithelial")][["stage", "bimodality_coefficient"]].set_index("stage")
    bc_q2_exit = q2[(q2["substrate"] == "he2022") & (q2["view"] == "exiting_tip")][["stage", "bimodality_coefficient"]].set_index("stage")
    bc_q21_all = out[out["view"] == "all_in_pool"][["stage", "bimodality_coefficient"]].set_index("stage")
    bc_q21_exit = out[out["view"] == "exiting_tip_in_pool"][["stage", "bimodality_coefficient"]].set_index("stage")
    cmp = pd.DataFrame({
        "Q2 all_epithelial": bc_q2_all["bimodality_coefficient"],
        "Q2 exiting_tip":    bc_q2_exit["bimodality_coefficient"],
        "Q2.1 all_in_pool":  bc_q21_all["bimodality_coefficient"],
        "Q2.1 exit_in_pool": bc_q21_exit["bimodality_coefficient"],
    }).reindex(HE_STAGE_ORDER).round(4)
    print(cmp.to_string())
    print(f"\n(BC threshold for switch-like: {BC_THRESHOLD:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
