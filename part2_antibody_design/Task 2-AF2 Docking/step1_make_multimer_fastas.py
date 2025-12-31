#!/usr/bin/env python3
"""
AbDesign: Multimer FASTA Generator
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script prepares input files for AlphaFold-Multimer. It combines candidate 
VH/VL sequences with a fixed Antigen sequence into individual FASTA files 
formatted for complex prediction (chain breaks indicated by separate entries).

Features:
---------
- Combines Antibody (VH+VL) and Antigen sequences
- Generates compliant FASTA headers for ColabFold
- Batch processing from CSV input

Requirements:
-------------
- pandas

Usage (example: tslp antigen):
------
    python step1_make_multimer_fastas.py --candidates_csv candidates.csv \
      --antigen_fasta tslp.fasta --out_dir fastas_multimer

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import os
import pandas as pd

def read_fasta_one(path: str) -> str:
    """Read a single sequence from a FASTA file."""
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).replace(" ", "").upper()

def main():
    parser = argparse.ArgumentParser(
        description="Generate Antibody:Antigen complex FASTAs for AlphaFold-Multimer."
    )
    parser.add_argument("--candidates_csv", required=True, help="CSV with columns: id, vh_seq, vl_seq")
    parser.add_argument("--antigen_fasta", required=True, help="FASTA file containing the antigen sequence")
    parser.add_argument("--out_dir", default="fastas_multimer", help="Output directory for FASTA files")
    args = parser.parse_args()

    # 1. Setup
    os.makedirs(args.out_dir, exist_ok=True)
    antigen_seq = read_fasta_one(args.antigen_fasta)
    
    print(f"[INFO] Antigen length: {len(antigen_seq)} residues")

    # 2. Load Candidates
    df = pd.read_csv(args.candidates_csv)
    required_cols = {"id", "vh_seq", "vl_seq"}
    missing = required_cols - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] Missing columns in CSV: {missing}")

    # 3. Generate Files
    count = 0
    for _, r in df.iterrows():
        cid = str(r["id"])
        vh = str(r["vh_seq"]).strip().upper()
        vl = str(r["vl_seq"]).strip().upper()

        out_fa = os.path.join(args.out_dir, f"{cid}.fasta")
        
        # Format: 3 separate chains (VH, VL, Antigen)
        # ColabFold handles multiple > entries as separate chains in multimer mode
        with open(out_fa, "w") as f:
            f.write(f">{cid}|VH\n{vh}\n")
            f.write(f">{cid}|VL\n{vl}\n")
            f.write(f">{cid}|Ag\n{antigen_seq}\n")
        count += 1

    print(f"[SUCCESS] Wrote {count} multimer FASTAs to: {args.out_dir}")

if __name__ == "__main__":
    main()