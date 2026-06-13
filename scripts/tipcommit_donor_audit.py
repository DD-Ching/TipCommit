#!/usr/bin/env python3
"""tipcommit donor-audit - quantify per-stage donor confounding in He 2022
and Cao 2020.

Background: q1_summary acknowledged "with 10 donors across 6 stages, some
stages have 1-2 donors" but did not quantify how donor-confounded the
per-stage trajectory actually is. This pass reads the per-cell scores
CSV produced by tipcommit_q2_mvp.py and emits, per (substrate, stage):

  - n_donors with >=1 cell at that stage
  - dominant donor's share of stage cells
  - n_donors after a 5%-of-stage minimum threshold
  - whether the candidate-pool subset is single-donor

Pure pandas. No Census re-fetch. No new dependencies.

Output:
  metadata/donor_audit_per_stage.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
META = REPO_ROOT / "metadata"
CANDIDATE_POOL_HE = "epithelial cell of lung"
CANDIDATE_POOL_CAO = "epithelial cell of lower respiratory tract"


def audit_substrate(df: pd.DataFrame, substrate: str, pool_ct: str,
                    stage_order: list[str]) -> pd.DataFrame:
    sub = df[df["substrate"] == substrate].copy()
    rows = []
    for stage in stage_order:
        st = sub[sub["development_stage"].astype(str) == stage]
        if len(st) == 0:
            continue
        donor_counts = st["donor_id"].value_counts()
        n_donors_any = (donor_counts > 0).sum()
        threshold_5pct = max(1, int(0.05 * len(st)))
        n_donors_5pct = (donor_counts >= threshold_5pct).sum()
        dominant_share = donor_counts.iloc[0] / len(st) if len(donor_counts) else float("nan")

        pool = st[st["cell_type"] == pool_ct]
        pool_donors = pool["donor_id"].value_counts()
        n_pool_donors_any = (pool_donors > 0).sum()
        pool_threshold = max(1, int(0.05 * len(pool))) if len(pool) else 1
        n_pool_donors_5pct = (pool_donors >= pool_threshold).sum() if len(pool) else 0

        rows.append({
            "substrate": substrate,
            "stage": stage,
            "n_total_epithelial": len(st),
            "n_pool": len(pool),
            "n_donors_any_cell": int(n_donors_any),
            "n_donors_geq5pct_share": int(n_donors_5pct),
            "dominant_donor_share": round(float(dominant_share), 3),
            "is_single_donor_dominant": bool(dominant_share > 0.95) if dominant_share == dominant_share else False,
            "n_pool_donors_any_cell": int(n_pool_donors_any),
            "n_pool_donors_geq5pct_share": int(n_pool_donors_5pct),
        })
    return pd.DataFrame(rows)


def main() -> int:
    in_path = META / "q2_per_cell_scores.csv"
    if not in_path.exists():
        print(f"missing input: {in_path}", file=sys.stderr)
        return 1
    df = pd.read_csv(in_path)
    print(f"read {in_path.relative_to(REPO_ROOT)}: {len(df):,} cells")

    he_stages = sorted(
        df[df["substrate"] == "he2022"]["development_stage"].astype(str).unique(),
        key=lambda s: int(s.split()[0].rstrip("thrdsnt")) if s[0].isdigit() else 999,
    )
    cao_stages = sorted(
        df[df["substrate"] == "cao2020"]["development_stage"].astype(str).unique(),
        key=lambda s: int(s.split()[0].rstrip("thrdsnt")) if s[0].isdigit() else 999,
    )

    he = audit_substrate(df, "he2022", CANDIDATE_POOL_HE, he_stages)
    cao = audit_substrate(df, "cao2020", CANDIDATE_POOL_CAO, cao_stages)
    out = pd.concat([he, cao], axis=0).reset_index(drop=True)

    out_path = META / "donor_audit_per_stage.csv"
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path.relative_to(REPO_ROOT)} ({out_path.stat().st_size} bytes)")

    print("\n=== Donor audit per (substrate, stage) ===")
    print(out.to_string(index=False))

    print("\n=== Summary ===")
    he_single = out[(out["substrate"] == "he2022")
                    & (out["dominant_donor_share"] > 0.95)]
    print(f"He 2022 stages with >95% from a single donor: "
          f"{len(he_single)} of {len(out[out['substrate']=='he2022'])} "
          f"({', '.join(he_single['stage'].tolist())})")
    he_balanced = out[(out["substrate"] == "he2022")
                      & (out["n_donors_geq5pct_share"] >= 2)]
    print(f"He 2022 stages with >=2 donors at >=5% share: "
          f"{len(he_balanced)} of {len(out[out['substrate']=='he2022'])} "
          f"({', '.join(he_balanced['stage'].tolist())})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
