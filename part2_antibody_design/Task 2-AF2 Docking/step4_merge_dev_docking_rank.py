#!/usr/bin/env python3
"""
AbDesign: Final Ranking (Developability + Docking)
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script merges the Developability Scores (CDS) calculated in Part II 
with the Docking Metrics (iPTM, iPAE) from step3_parse_af2_and_interface.py It computes a 
Final Score to rank candidates for experimental validation.

Features:
---------
- Merges disparate CSV datasets by Candidate ID
- Normalizes metrics (MinMax scaling) to enable linear combination
- Computes Final Score = w_dev*CDS + w_iptm*iPTM - w_pae*iPAE
- Exports final rankings with provenance

Requirements:
-------------
- pandas
- numpy

Usage:
------
    python step4_merge_dev_docking_rank.py \
      --dev_csv candidates_ranked.csv \
      --dock_csv summary_docking.csv \
      --out_csv final_selection.csv

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import pandas as pd

def norm_minmax(series: pd.Series) -> pd.Series:
    """Normalize a series to [0,1]."""
    x = pd.to_numeric(series, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series([0.0]*len(series), index=series.index)
    lo, hi = float(np.nanmin(x.values)), float(np.nanmax(x.values))
    if hi == lo: return pd.Series([0.0]*len(series), index=series.index)
    return (x - lo) / (hi - lo)

def main():
    parser = argparse.ArgumentParser(description="Merge scores and rank candidates.")
    parser.add_argument("--dev_csv", required=True, help="Part II output (CDS scores)")
    parser.add_argument("--dock_csv", required=True, help="Part III output (Docking metrics)")
    parser.add_argument("--out_csv", required=True, help="Final ranking CSV")
    
    parser.add_argument("--dev_col", default="DCS", help="Column name for Developability Score")
    parser.add_argument("--min_iptm", type=float, default=0.2, help="Filter: minimum docking iPTM")
    
    # Weights for Final Score
    parser.add_argument("--w_dev", type=float, default=1.0, help="Weight: Developability")
    parser.add_argument("--w_iptm", type=float, default=1.5, help="Weight: Docking iPTM")
    parser.add_argument("--w_pae", type=float, default=0.5, help="Weight: Interface PAE (penalty)")
    
    args = parser.parse_args()

    # 1. Load Data
    dev = pd.read_csv(args.dev_csv)
    dock = pd.read_csv(args.dock_csv)

    # 2. ID Harmonization
    # Ensure join keys match. Part II usually exports 'id', Part III 'candidate_id'
    if "candidate_id" not in dev.columns and "id" in dev.columns:
        dev = dev.rename(columns={"id": "candidate_id"})
    
    # Clean IDs (remove extra suffices like _complex or _fv if present to ensure match)
    # This assumes IDs in both files share a common base (e.g. tezepelumab_var_0001)
    # For now, we assume simple exact matching or simple substring containment.
    
    # 3. Merge
    print(f"[INFO] Merging {len(dev)} dev records with {len(dock)} docking records...")
    # Attempt cleanup if standard merge fails
    if not set(dev["candidate_id"]) & set(dock["candidate_id"]):
         print("[WARN] ID mismatch detected. Attempting to strip suffixes...")
         dev["candidate_id"] = dev["candidate_id"].astype(str).str.replace("_fv", "").str.replace("_complex", "")
         dock["candidate_id"] = dock["candidate_id"].astype(str).str.replace("_fv", "").str.replace("_complex", "")

    m = dev.merge(dock, on="candidate_id", how="inner")
    
    if len(m) == 0:
        raise SystemExit("[ERROR] Merge produced 0 rows. Check candidate IDs.")

    # 4. Filter
    if args.min_iptm:
        m = m[m["iptm"] >= args.min_iptm].copy()
        print(f"[INFO] {len(m)} candidates remain after iPTM >= {args.min_iptm} filter.")

    # 5. Normalize Scores
    m["dev_norm"] = norm_minmax(m[args.dev_col])
    m["iptm_norm"] = norm_minmax(m["iptm"])
    # PAE is 'lower is better', so we subtract it. Normalize normally first.
    m["pae_norm"] = norm_minmax(m["mean_interface_pae"])

    # 6. Final Composite Score
    # Score = (Dev * w1) + (iPTM * w2) - (PAE * w3)
    m["final_score"] = (
        args.w_dev * m["dev_norm"] + 
        args.w_iptm * m["iptm_norm"] - 
        args.w_pae * m["pae_norm"]
    )

    # 7. Rank & Export
    m = m.sort_values("final_score", ascending=False)
    m["final_rank"] = range(1, len(m) + 1)

    # Select columns
    cols = ["final_rank", "candidate_id", "final_score", 
            args.dev_col, "iptm", "mean_interface_pae", 
            "plddt", "pdb_path"]
    
    # Add optional columns if present
    for c in ["liability_risk_cdr", "solubility_score"]:
        if c in m.columns: cols.append(c)

    m[cols].to_csv(args.out_csv, index=False)
    
    print(f"[SUCCESS] Final ranking saved to {args.out_csv}")
    print("[INFO] Top 5 Candidates:")
    print(m[cols].head(5).to_string(index=False))

if __name__ == "__main__":
    main()