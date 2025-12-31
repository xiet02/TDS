#!/usr/bin/env python3
"""
AbDesign: CDR Liability Scorer
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script identifies chemical liability motifs specifically within the CDRs 
of scFv sequences. It splits scFvs into VH/VL domains and applies weighted 
penalties for motifs such as glycosylation sites, deamidation, and oxidation.

Features:
---------
- Automatic VH/VL splitting via linker detection
- Identification of N-glycan, Deamidation, Isomerization, and Oxidation motifs
- Weighted Liability Risk calculation [cdr only]
- Robust handling of varied linker lengths

Requirements:
-------------
- pandas
- biopython

Usage:
------
    python step3_liability_scoring_cdr.py input.fasta --out liabilities.csv

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import math
import re
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# ============================================================================
# Configuration & Weights
# ============================================================================

DG_WEIGHT = 1.5
TAU = 12.0

WEIGHTS = {
    "cdr_nglyco_NXS_T": 3.0,
    "cdr_isomer_DG":    DG_WEIGHT,
    "cdr_deamid_NG":    1.0,
    "cdr_deamid_NS":    1.0,
    "cdr_deamid_NT":    1.0,
    "cdr_deamid_NN":    1.0,
    "cdr_cleavage_DP":  1.5,
    "cdr_oxid_M":       0.5,
    "cdr_oxid_W":       0.5,
}

# CDR Definitions (0-based indices)
# example: teztuzumab (Kabat numbering)
L1 = list(range(23 - 1, 33))
L2 = list(range(49 - 1, 55))
L3 = list(range(88 - 1, 98))

H1 = list(range(26 - 1, 35))
H2 = list(range(50 - 1, 66))
H3 = list(range(98 - 1, 111))

CDR_POS_H = sorted(set(H1 + H2 + H3))
CDR_POS_L = sorted(set(L1 + L2 + L3))

# ============================================================================
# Motif Detection Logic
# ============================================================================

def count_motifs(seq: str) -> dict:
    """Scan sequence for liability motifs."""
    c = defaultdict(int)

    # N-glycosylation NXS/T (X != P)
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in ("S", "T"):
            c["cdr_nglyco_NXS_T"] += 1

    c["cdr_isomer_DG"]   += len(re.findall("DG", seq))
    c["cdr_deamid_NG"]   += len(re.findall("NG", seq))
    c["cdr_deamid_NS"]   += len(re.findall("NS", seq))
    c["cdr_deamid_NT"]   += len(re.findall("NT", seq))
    c["cdr_deamid_NN"]   += len(re.findall("NN", seq))
    c["cdr_cleavage_DP"] += len(re.findall("DP", seq))
    c["cdr_oxid_M"]      += seq.count("M")
    c["cdr_oxid_W"]      += seq.count("W")

    return c

def compute_liability_risk_cdr(counts: dict, length_cdr: int) -> float:
    """Calculate normalized risk score [0,1] based on motif density."""
    weighted = 0.0
    for k, w in WEIGHTS.items():
        weighted += w * counts.get(k, 0)

    L = max(length_cdr, 1)
    density = 100.0 * weighted / L
    risk = 1.0 - math.exp(-density / TAU)
    return round(max(0.0, min(1.0, risk)), 6)

# ============================================================================
# Sequence Splitting
# ============================================================================

def split_scfv_vh_vl(seq: str, min_linker_len: int = 12):
    """
    Split scFv into VH and VL by finding the GS-linker.
    
    Args:
        seq (str): scFv sequence
        min_linker_len (int): Min length for linker detection
        
    Returns:
        tuple: (vh, linker, vl)
    """
    # Find longest run of G/S/P
    best = (None, None, 0) # start, end, length
    i, n = 0, len(seq)
    allowed = set("GSP")
    
    while i < n:
        if seq[i] in allowed:
            j = i
            while j < n and seq[j] in allowed:
                j += 1
            if (j - i) > best[2]:
                best = (i, j, j - i)
            i = j
        else:
            i += 1

    if best[2] < min_linker_len:
        raise ValueError("Linker too short/not found")

    s, e = best[0], best[1]
    vh = seq[:s]
    linker = seq[s:e]
    vl = seq[e:]
    
    return vh, linker, vl

def extract_cdr_by_positions(domain_seq: str, positions: list) -> str:
    return "".join([domain_seq[i] for i in positions if 0 <= i < len(domain_seq)])

# ============================================================================
# Main Entry Point
# ============================================================================

def main(fasta, out_csv, min_linker_len):
    rows = []
    print(f"[INFO] Scanning liabilities in {fasta}...")

    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper().replace(" ", "").replace("\n", "")
        if not seq: continue

        row = {"id": rec.id}

        try:
            # 1. Split scFv
            vh, linker, vl = split_scfv_vh_vl(seq, min_linker_len)
            
            # 2. Extract CDRs
            vh_cdr = extract_cdr_by_positions(vh, CDR_POS_H)
            vl_cdr = extract_cdr_by_positions(vl, CDR_POS_L)
            
            # 3. Count Motifs
            c_h = count_motifs(vh_cdr)
            c_l = count_motifs(vl_cdr)
            
            # Merge counts
            counts = defaultdict(int)
            for d in [c_h, c_l]:
                for k,v in d.items(): counts[k] += v
            
            length_cdr = len(vh_cdr) + len(vl_cdr)
            risk = compute_liability_risk_cdr(counts, length_cdr)
            
            row.update({
                "split_ok": 1,
                "length_cdr": length_cdr,
                "liability_risk_cdr": risk
            })
            row.update(counts)

        except Exception:
            # Fallback for failed splits
            c = count_motifs(seq)
            row.update({
                "split_ok": 0,
                "length_cdr": len(seq),
                "liability_risk_cdr": compute_liability_risk_cdr(c, len(seq))
            })
            row.update(c)
        
        rows.append(row)

    # Save Results
    df = pd.DataFrame(rows).fillna(0)
    
    # Ensure columns exist
    for k in WEIGHTS.keys():
        if k not in df.columns: df[k] = 0
            
    df.to_csv(out_csv, index=False)
    print(f"[SUCCESS] Wrote {out_csv} ({len(df)} sequences)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="CDR-only antibody liability scoring (Kabat-indexed)."
    )
    ap.add_argument("fasta", help="Input FASTA (scFv: VH-linker-VL)")
    ap.add_argument("-o", "--out", default="liabilities_cdr_v2.csv")
    ap.add_argument("--min_linker_len", type=int, default=12, help="Min linker length")
    args = ap.parse_args()

    main(args.fasta, args.out, args.min_linker_len)