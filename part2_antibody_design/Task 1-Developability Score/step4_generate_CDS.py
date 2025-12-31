#!/usr/bin/env python3
"""
AbDesign: Composite Developability Score (CDS) Generator
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script aggregates data from structure prediction (AF2), solubility proxies,
and liability scoring to compute the final Composite Developability Score (CDS).
It ranks candidates and exports the top performers for experimental validation.

Features:
---------
- Data merging from multi-stage pipeline outputs
- ID normalization and error checking
- Hard filtering of non-viable candidates (low pLDDT, glycosylation)
- Weighted CDS calculation [cite: 56-58]
- Diversity-aware selection of top candidates

Requirements:
-------------
- pandas
- biopython

Usage:
------
    python step4_generate_CDS.py

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import re
import pandas as pd
from Bio import SeqIO

# ============================================================================
# Configuration
# ============================================================================

# Input Files, modify paths/filenames as needed
AF2_FILE = "summary/af2_summary_rank001.csv"
SOL_FILE = "solubility_proxy.csv"
LIA_FILE = "liabilities_cdr_v2.csv"
FASTA_FILE = "data/tezepelumab_lib_fv2.fasta"

# Output Files, modify paths/filenames as needed
OUT_ALL = "candidates_ranked.csv"
OUT_TOP = "top20.fasta"

# ============================================================================
# ID Normalization
# ============================================================================

def canon_id(x: str) -> str:
    """Standardize variant IDs across different pipeline steps."""
    if pd.isna(x): return None
    s = str(x)
    
    # Matches: tezepelumab_var_XXXX_fv
    m = re.search(r"(tezepelumab_var_\d{4}_fv)", s)
    if m: return m.group(1)

    # Matches: tezepelumab_var_XXXX
    m = re.search(r"(tezepelumab_var_\d{4})", s)
    if m: return m.group(1) + "_fv"

    return s

# ============================================================================
# Main Logic
# ============================================================================

def main():
    print("[INFO] Loading datasets...")
    
    # 1. Load Data
    try:
        af2 = pd.read_csv(AF2_FILE)
        sol = pd.read_csv(SOL_FILE)
        lia = pd.read_csv(LIA_FILE)
    except FileNotFoundError as e:
        print(f"[ERROR] Missing input file: {e}")
        return

    # 2. Normalize IDs
    af2["id"] = af2["id"].map(canon_id)
    sol["id"] = sol["id"].map(canon_id)
    lia["id"] = lia["id"].map(canon_id)

    # 3. Merge Data
    df = af2.merge(sol, on="id", how="inner").merge(lia, on="id", how="inner")
    print(f"[INFO] Merged {len(df)} variants common to all datasets.")

    if len(df) == 0:
        raise SystemExit("[ERROR] Merge resulted in 0 rows. Check ID formatting.")

    # 4. Hard Filters
    # pLDDT > 80, Solubility > 0.45, No Glycosylation
    df = df[
        (df["mean_plddt"] >= 80) &
        (df["fw_plddt"] >= 88) &
        (df["solubility_score"] >= 0.45) &
        (df["cdr_nglyco_NXS_T"] == 0)
    ].copy()
    print(f"[INFO] Variants passing hard filters: {len(df)}")

    # 5. Calculate Scores
    # Normalize pLDDT to 0-1
    df["mean_plddt_n"] = df["mean_plddt"] / 100.0
    df["fw_plddt_n"]   = df["fw_plddt"] / 100.0
    df["cdr_plddt_n"]  = df["cdr_plddt"] / 100.0

    # Structural Score
    df["struct_score"] = 0.5*df["mean_plddt_n"] + 0.3*df["cdr_plddt_n"] + 0.2*df["fw_plddt_n"]
    
    # Stability Score Proxy
    df["stability_score"] = df["fw_plddt_n"]

    # Composite Developability Score (CDS/DCS)
    # Weights: 30% Structure, 25% Solubility, 25% Low Liability, 20% Stability
    df["DCS"] = (
        0.30 * df["struct_score"] +
        0.25 * df["solubility_score"] +
        0.25 * (1.0 - df["liability_risk_cdr"]) +
        0.20 * df["stability_score"]
    ) * 100.0

    # Sort and Save
    df = df.sort_values("DCS", ascending=False)
    df.to_csv(OUT_ALL, index=False)
    print(f"[SUCCESS] Ranked candidates saved to {OUT_ALL}")

    # 6. Select Top Diverse Candidates
    # Bucket by features to ensure diversity in top selection
    df["bucket"] = (
        pd.cut(df["solubility_score"], bins=5, labels=False).astype(str) + "_" +
        pd.cut(df["liability_risk_cdr"], bins=5, labels=False).astype(str) + "_" +
        pd.cut(df["cdr_plddt"], bins=5, labels=False).astype(str)
    )

    # Pick top 2 from each bucket, up to 20 total
    top = df.groupby("bucket", as_index=False, sort=False).head(2).head(20)

    # 7. Export Top FASTA
    seqs = {rec.id: rec for rec in SeqIO.parse(FASTA_FILE, "fasta")}
    count = 0
    with open(OUT_TOP, "w") as f:
        for _, r in top.iterrows():
            rec = seqs.get(r["id"])
            if rec:
                SeqIO.write(rec, f, "fasta")
                count += 1

    print(f"[SUCCESS] Exported {count} top candidates to {OUT_TOP}")

if __name__ == "__main__":
    main()