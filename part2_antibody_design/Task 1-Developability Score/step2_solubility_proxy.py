#!/usr/bin/env python3
"""
AbDesign: Solubility Proxy Calculator
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script calculates sequence-based solubility proxies for antibody candidates.
It implements the scoring logic described in papers, penalizing 
high hydrophobicity, aromatic density, and extreme charge deviations.

Features:
---------
- Kyte-Doolittle hydrophobicity analysis
- Net charge calculation (pH 7.4 proxy)
- Aromatic and hydrophobic fraction assessment
- Normalized solubility score (0.0 - 1.0)

Requirements:
-------------
- pandas
- numpy
- biopython

Usage:
------
    python step2_solubility_proxy.py library.fasta --out solubility_scores.csv

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import math
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

# ============================================================================
# Constants & Scales
# ============================================================================

# Kyte-Doolittle hydrophobicity scale
KD: Dict[str, float] = {
    "A": 1.8,  "C": 2.5,  "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5,  "K": -3.9, "L": 3.8,
    "M": 1.9,  "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2,  "W": -0.9, "Y": -1.3
}

POS = set("KRH")
NEG = set("DE")
AROM = set("FWY")
HYDROPHOBIC_AA = set("AVILMFWYC")

# ============================================================================
# Scoring Logic
# ============================================================================

def clean_seq(seq: str) -> str:
    """Normalize sequence to uppercase standard amino acids."""
    seq = seq.strip().upper()
    allowed = set(KD.keys())
    return "".join([aa for aa in seq if aa in allowed])

def mean_hydrophobicity(seq: str) -> float:
    vals = [KD[aa] for aa in seq]
    return float(np.mean(vals)) if vals else float("nan")

def net_charge_proxy(seq: str) -> float:
    """Approximate net charge at pH 7.4 (Histidine = +0.1)."""
    k = seq.count("K")
    r = seq.count("R")
    h = seq.count("H")
    d = seq.count("D")
    e = seq.count("E")
    return (k + r + 0.1 * h) - (d + e)

def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))

def solubility_score(seq: str) -> Tuple[float, Dict[str, float]]:
    """
    Calculate solubility score [0,1].
    
    Logic:
    1. Calculate phys-chem properties.
    2. Normalize properties into penalty terms (0=good, 1=bad).
    3. Weighted sum of penalties.
    4. Final Score = 1 - Total Penalty.
    """
    L = len(seq)
    mh = mean_hydrophobicity(seq)
    nc = net_charge_proxy(seq)
    
    fa = sum(1 for aa in seq if aa in AROM) / L if L > 0 else 0
    fh = sum(1 for aa in seq if aa in HYDROPHOBIC_AA) / L if L > 0 else 0

    # 1. Hydrophobicity Penalty: mh > 0.3 is bad
    p_hydro = sigmoid((mh - 0.3) * 2.0)

    # 2. Aromatic Penalty: > 8% is sticky
    p_arom = sigmoid((fa - 0.08) * 25.0)

    # 3. Charge Penalty: Extreme charge density is bad
    nc_norm = (nc / max(L, 1)) * 100.0
    p_charge = sigmoid((abs(nc_norm) - 8.0) * 0.8)

    # 4. Hydrophobic Fraction Penalty: > 45% is aggregation prone
    p_hfrac = sigmoid((fh - 0.45) * 10.0)

    # Composite Penalty (Weights derived from historical data)
    penalty = 0.45 * p_hydro + 0.25 * p_arom + 0.20 * p_charge + 0.10 * p_hfrac

    score = max(0.0, min(1.0, 1.0 - penalty))

    feats = {
        "length": L,
        "mean_hydrophobicity": mh,
        "net_charge_proxy": nc,
        "frac_aromatic": fa,
        "solubility_score": score,
    }
    return score, feats

# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate solubility proxy scores for antibody sequences."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("-o", "--out", default="solubility_proxy.csv", help="Output CSV path")
    args = parser.parse_args()

    rows: List[Dict[str, float]] = []
    
    print(f"[INFO] Scoring sequences from {args.fasta}...")
    for rec in SeqIO.parse(args.fasta, "fasta"):
        seq = clean_seq(str(rec.seq))
        if not seq:
            continue
        _, feats = solubility_score(seq)
        feats["id"] = rec.id
        rows.append(feats)

    df = pd.DataFrame(rows)
    df = df.sort_values("solubility_score", ascending=False)
    df.to_csv(args.out, index=False)
    
    print(f"[SUCCESS] Wrote {args.out} with {len(df)} sequences.")

if __name__ == "__main__":
    main()