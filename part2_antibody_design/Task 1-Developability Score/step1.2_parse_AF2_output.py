#!/usr/bin/env python3
"""
AbDesign: AlphaFold2 Output Parser
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script performs detailed parsing of AlphaFold2 output JSONs. unlike the
simple summary in step 1.1, this script calculates domain-specific confidence
metrics (Framework vs CDR pLDDT) required for the Composite Developability Score.

Features:
---------
- Extracts global confidence metrics (pTM, iPTM)
- Calculates region-specific pLDDT (Framework, CDR, Min/Max)
- Aggregates results into a Pandas DataFrame
- Supports ranking filtering

Requirements:
-------------
- pandas
- numpy

Usage:
------
    python step1.2_parse_AF2_output.py af2_results/ --rank 1

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import glob
import json
import os
import numpy as np
import pandas as pd

# ============================================================================
# Metric Calculation
# ============================================================================

def job_from_filename(path: str) -> str:
    """Extract clean job ID from filename."""
    base = os.path.basename(path)
    if "_scores_rank_001_" in base:
        return base.split("_scores_rank_001_")[0]
    if "_scores_" in base:
        return base.split("_scores_")[0]
    return os.path.splitext(base)[0]

def compute_plddt_metrics(plddt_list):
    """
    Calculate regional pLDDT scores.
    
    Assumes standard scFv length and approximates regions:
    - Framework: N-term and C-term ends
    - CDRs: Central region (rough approximation for sorting)
    """
    p = np.array(plddt_list, dtype=float)
    mean_plddt = float(p.mean())
    min_plddt  = float(p.min())
    max_plddt  = float(p.max())

    # Framework proxy: Ends of sequence (first and last 80 residues)
    # If sequence is short, uses whole sequence
    fw = np.concatenate([p[:80], p[-80:]]) if len(p) >= 160 else p
    fw_plddt = float(fw.mean())

    # CDR proxy: Central band (approx 30%-70% of seq)
    mid_start = max(0, int(len(p) * 0.30))
    mid_end   = min(len(p), int(len(p) * 0.70))
    cdr = p[mid_start:mid_end] if mid_end > mid_start else p
    cdr_plddt = float(cdr.mean())

    return mean_plddt, fw_plddt, cdr_plddt, min_plddt, max_plddt

# ============================================================================
# Main Parsing Logic
# ============================================================================

def main(result_dir: str, out_csv: str, rank: int):
    pattern = os.path.join(result_dir, f"*_scores_rank_{rank:03d}_*.json")
    files = sorted(glob.glob(pattern))
    
    if not files:
        raise SystemExit(f"[ERROR] No score json files matched: {pattern}")

    rows = []
    print(f"[INFO] Processing {len(files)} score files from {result_dir}...")

    for fpath in files:
        with open(fpath, "r") as f:
            d = json.load(f)

        job = job_from_filename(fpath)
        plddt = d.get("plddt", None)
        
        if plddt is None:
            continue

        mean_plddt, fw_plddt, cdr_plddt, min_plddt, max_plddt = compute_plddt_metrics(plddt)

        rows.append({
            "id": job,
            "score_file": os.path.basename(fpath),
            "mean_plddt": mean_plddt,
            "fw_plddt": fw_plddt,
            "cdr_plddt": cdr_plddt,
            "min_plddt": min_plddt,
            "max_plddt": max_plddt,
            "ptm": d.get("ptm", None),
            "iptm": d.get("iptm", None),
            "iptm_plus_ptm": d.get("iptm+ptm", None),
            "ranking_confidence": d.get("ranking_confidence", None),
            "max_pae": d.get("max_pae", None),
        })

    df = pd.DataFrame(rows)
    # Deduplicate and sort
    df = df.drop_duplicates(subset=["id"]).sort_values("mean_plddt", ascending=False)
    
    # Save
    df.to_csv(out_csv, index=False)
    print(f"[SUCCESS] Wrote extended summary: {out_csv} ({len(df)} sequences)")
    

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Generate detailed AF2 metrics CSV from ColabFold output."
    )
    ap.add_argument("result_dir", help="Folder containing AF2 JSON outputs")
    ap.add_argument("-o", "--out", default=None, help="Output CSV path (default: <result_dir>/af2_summary_rankXXX.csv)")
    ap.add_argument("--rank", type=int, default=1, help="Rank index to process (default: 1)")
    
    args = ap.parse_args()

    out_csv = args.out or os.path.join(args.result_dir, f"af2_summary_rank{args.rank:03d}.csv")
    main(args.result_dir, out_csv, args.rank)