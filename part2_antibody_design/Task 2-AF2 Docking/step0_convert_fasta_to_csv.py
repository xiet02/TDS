#!/usr/bin/env python3
"""
AbDesign: Fv Sequence Parser & Splitter
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script processes a FASTA file containing scFv sequences (VH-Linker-VL).
It automatically detects common linker patterns to split the sequences into 
distinct Heavy (VH) and Light (VL) chain components required for multimer modeling.

Features:
---------
- Heuristic detection of GGGGS-type linkers
- Validation of domain lengths (>90 AA)
- CSV export for batch processing

Requirements:
-------------
- python >= 3.8

Usage:
------
    python step0_convert_fasta_to_csv.py --in_fasta library.fasta --out_csv candidates.csv

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import re
import csv
import argparse
from typing import List, Tuple

# Common flexible linker patterns seen in scFv constructs
LINKER_CANDIDATES = [
    "GGGGSGGGGSGGGGSS",
    "GGGGSGGGGSGGGGS",
    "GGGGS" * 3,            # Standard 15-mer
    "GGGGS" * 4,            # 20-mer
]

def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Read FASTA and return list of (header, sequence) tuples."""
    records = []
    cur_id = None
    cur_seq = []
    
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq).replace(" ", "").upper()))
                cur_id = line[1:].split()[0] # Take first word as ID
                cur_seq = []
            else:
                cur_seq.append(line)
    
    if cur_id is not None:
        records.append((cur_id, "".join(cur_seq).replace(" ", "").upper()))
    return records

def split_fv(seq: str) -> Tuple[str, str, str]:
    """
    Split scFv into VH and VL based on linker detection.
    Returns: (VH, VL, Linker)
    """
    # 1. Try exact known linker candidates
    for lk in LINKER_CANDIDATES:
        if lk in seq:
            vh, vl = seq.split(lk, 1)
            return vh, vl, lk

    # 2. Fallback: regex for (GGGGS){3,} plus optional tail S/SS
    # Captures standard flexible linkers with slight variations
    m = re.search(r"(?:GGGGS){3,}S{0,2}", seq)
    if m:
        lk = m.group(0)
        vh = seq[:m.start()]
        vl = seq[m.end():]
        return vh, vl, lk

    raise ValueError("No recognizable Fv linker found")

def main():
    parser = argparse.ArgumentParser(
        description="Split scFv FASTA into VH/VL columns for multimer docking."
    )
    parser.add_argument("--in_fasta", required=True, help="Input FASTA with one Fv per entry")
    parser.add_argument("--out_csv", default="candidates_split.csv", help="Output CSV path")
    args = parser.parse_args()

    records = read_fasta(args.in_fasta)
    rows = []

    print(f"[INFO] Processing {len(records)} sequences...")

    for rid, seq in records:
        try:
            vh, vl, lk = split_fv(seq)
            # Sanity check: Ab domains are usually ~110-120 AA
            if len(vh) < 90 or len(vl) < 90:
                print(f"[WARN] {rid}: Suspicious domain lengths (VH={len(vh)}, VL={len(vl)}). Skipping.")
                continue
            
            rows.append({"id": rid, "vh_seq": vh, "vl_seq": vl, "linker": lk})
        except ValueError as e:
            print(f"[WARN] {rid}: {e}")

    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "vh_seq", "vl_seq", "linker"])
        w.writeheader()
        w.writerows(rows)

    print(f"[SUCCESS] Wrote {len(rows)} valid candidates to {args.out_csv}")

if __name__ == "__main__":
    main()