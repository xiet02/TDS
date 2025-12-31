#!/usr/bin/env python3
"""
AbDesign: Docking Metric Parser (PAE & Interface)
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script parses the complex JSON and PDB outputs from AlphaFold-Multimer.
It calculates critical interface quality metrics, specifically the 
Interface Predicted Aligned Error (iPAE), which is key for ranking binders.

Features:
---------
- Extracts global metrics (pLDDT, ptm, iptm)
- Computes mean Interface PAE (lower is better) between Antibody and Antigen
- Auto-detects chain boundaries from PDB files
- Aggregates results into a summary CSV

Requirements:
-------------
- numpy

Usage:
------
    python step3_parse_af2_and_interface.py --results_dir results_docking --out_csv summary_docking.csv

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ============================================================================
# Helpers
# ============================================================================

def load_json(path: Path) -> Dict:
    with path.open("r") as f:
        return json.load(f)

def mean(x: List[float]) -> Optional[float]:
    return float(sum(x) / len(x)) if x else None

def pick_float_or_mean(obj, key: str) -> Optional[float]:
    """Extract scalar or mean of list from dictionary."""
    if key not in obj or obj[key] is None:
        return None
    v = obj[key]
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, list) and v and isinstance(v[0], (int, float)):
        return mean([float(z) for z in v])
    return None

def get_pae_matrix(obj: Dict) -> Optional[List[List[float]]]:
    """Robustly extract PAE matrix from various JSON schemas."""
    # Direct list
    for k in ["predicted_aligned_error", "pae"]:
        v = obj.get(k)
        if isinstance(v, list) and v and isinstance(v[0], list):
            return v
    # Nested dict (common in V1)
    if isinstance(obj.get("pae"), dict):
        v = obj["pae"].get("predicted_aligned_error")
        if isinstance(v, list): return v
    return None

def parse_chain_lengths_from_pdb(pdb_path: Path) -> List[Tuple[str, int]]:
    """Determine chain lengths from PDB ATOM records."""
    seen = set()
    order: List[str] = []
    counts: Dict[str, int] = {}

    with pdb_path.open("r") as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            if line[12:16].strip() != "CA": continue
            
            chain = (line[21].strip() or "_")
            resseq = line[22:26].strip()
            icode = (line[26].strip() or "")
            key = (chain, resseq, icode)
            
            if key in seen: continue
            seen.add(key)
            
            if chain not in counts:
                counts[chain] = 0
                order.append(chain)
            counts[chain] += 1

    return [(c, counts[c]) for c in order]

def mean_interface_pae_A_B(pae: List[List[float]], chain_lengths: List[Tuple[str, int]]) -> Optional[float]:
    """
    Calculate mean PAE for the interface (A <-> B interaction blocks).
    Assumes Chain A+B is Antibody (or A=Ab, B=Ag depending on complex).
    Here we calculate generic interface PAE between the first two distinct groups.
    """
    if not pae or not chain_lengths: return None

    # Determine blocks
    ranges = []
    start = 0
    for cid, n in chain_lengths:
        ranges.append((cid, start, start + n))
        start += n

    if len(ranges) < 2: return None

    # Heuristic: Interaction between Chain 0 (Ab part) and Chain 1 (Ag part)
    # Ideally, group VH+VL as "Ab" and Antigen as "Ag".
    # Since fasta generation wrote VH, VL, Ag, usually:
    # A=VH, B=VL, C=Ag.
    # We want interface between {A,B} and {C}.
    
    # Let's target the last chain (Antigen) vs the rest (Antibody)
    ag_range = ranges[-1]
    ab_ranges = ranges[:-1]
    
    vals = []
    
    # Ag residues
    _, c_s, c_e = ag_range
    
    for _, ab_s, ab_e in ab_ranges:
        # PAE[ab][ag]
        for i in range(ab_s, ab_e):
            vals.extend(pae[i][c_s:c_e])
        # PAE[ag][ab]
        for i in range(c_s, c_e):
            vals.extend(pae[i][ab_s:ab_e])

    return float(sum(vals) / len(vals)) if vals else None

def safe_float(x) -> Optional[float]:
    try:
        return float(x) if x is not None and not math.isnan(float(x)) else None
    except: return None

# ============================================================================
# Main Logic
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Extract AF2 docking metrics.")
    parser.add_argument("--results_dir", required=True, help="ColabFold output folder")
    parser.add_argument("--out_csv", required=True, help="Output summary CSV")
    parser.add_argument("--tag", default="complex", help="Filename tag to match (default: complex)")
    args = parser.parse_args()

    root = Path(args.results_dir)
    rows = []
    
    # 1. Group files by candidate
    # Pattern: candidate_id_complex_scores_rank_...
    # We scan recursively
    candidates = {}
    print(f"[INFO] Scanning {root}...")
    
    for p in root.rglob("*"):
        if not p.is_file(): continue
        # Heuristic ID extraction
        # Assumes filenames start with candidate ID
        if "_scores" in p.name:
            # teze_var_01_complex_scores_rank_001...
            cid = p.name.split("_scores")[0]
            candidates.setdefault(cid, {})["scores"] = candidates.get(cid, {}).get("scores", []) + [p]
        elif "_predicted_aligned_error" in p.name:
            cid = p.name.split("_predicted")[0]
            candidates.setdefault(cid, {})["pae"] = candidates.get(cid, {}).get("pae", []) + [p]
        elif p.suffix == ".pdb" and "rank_" in p.name:
             cid = p.name.split("_unrelaxed")[0].split("_relaxed")[0]
             candidates.setdefault(cid, {})["pdb"] = candidates.get(cid, {}).get("pdb", []) + [p]

    print(f"[INFO] Found {len(candidates)} candidates.")

    # 2. Process Best Rank
    for cid, files in sorted(candidates.items()):
        # Find Rank 001 files
        scores_f = sorted([f for f in files.get("scores", []) if "rank_001" in f.name])
        pae_f    = sorted([f for f in files.get("pae", []) if "rank_001" in f.name])
        pdb_f    = sorted([f for f in files.get("pdb", []) if "rank_001" in f.name])

        if not scores_f: continue
        
        # Load Data
        sj = load_json(scores_f[0])
        
        # Metrics
        plddt = pick_float_or_mean(sj, "plddt")
        ptm   = pick_float_or_mean(sj, "ptm")
        iptm  = pick_float_or_mean(sj, "iptm")

        # PAE Matrix
        pae_matrix = get_pae_matrix(sj)
        if pae_matrix is None and pae_f:
             pae_matrix = get_pae_matrix(load_json(pae_f[0]))

        # Interface PAE
        iface_pae = None
        chain_lens = []
        if pdb_f and pdb_f[0].exists():
            chain_lens = parse_chain_lengths_from_pdb(pdb_f[0])
            if pae_matrix:
                iface_pae = mean_interface_pae_A_B(pae_matrix, chain_lens)

        rows.append({
            "candidate_id": cid,
            "plddt": safe_float(plddt),
            "ptm": safe_float(ptm),
            "iptm": safe_float(iptm),
            "mean_interface_pae": safe_float(iface_pae),
            "chain_lengths": ";".join([f"{c}:{n}" for c,n in chain_lens]),
            "pdb_path": str(pdb_f[0]) if pdb_f else ""
        })

    # 3. Export
    keys = ["candidate_id", "plddt", "ptm", "iptm", "mean_interface_pae", "chain_lengths", "pdb_path"]
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print(f"[SUCCESS] Wrote metrics for {len(rows)} models to {args.out_csv}")

if __name__ == "__main__":
    main()