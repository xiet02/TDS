#!/usr/bin/env python3
"""
AbDesign: ColabFold Batch Executor
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script automates the execution of `colabfold_batch` for antibody-antigen 
docking. It iterates through generated FASTA files and runs structure prediction
using AlphaFold-Multimer.

Features:
---------
- Batch processing of directory inputs
- Configurable model parameters (recycles, seeds, models)
- Organized output structure (one folder per candidate)
- Dry-run mode for verification

Requirements:
-------------
- colabfold_batch (in system PATH)

Usage:
------
    python step2_run_colabfold_loop.py --fasta_dir fastas_multimer --results_dir results

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import glob
import os
import subprocess
from pathlib import Path

def run_command(cmd: list, dry_run: bool = False):
    """Print and execute a shell command."""
    print(" ".join(cmd))
    if dry_run:
        return
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run AlphaFold-Multimer on batch FASTAs.")
    parser.add_argument("--fasta_dir", default="fastas_multimer", help="Directory containing input FASTAs")
    parser.add_argument("--results_dir", default="results_docking", help="Directory for output results")
    parser.add_argument("--model_type", default="alphafold2_multimer_v3", help="ColabFold model type")
    parser.add_argument("--num_models", type=int, default=5, help="Number of models to generate")
    parser.add_argument("--num_recycle", type=int, default=3, help="Number of recycles (default: 3)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--gpu", action="store_true", help="Pass-through flag (implicit in ColabFold)")
    parser.add_argument("--dry_run", action="store_true", help="Print commands without executing")
    args = parser.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)

    # Find FASTAs
    fasta_paths = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fasta")))
    if not fasta_paths:
        raise SystemExit(f"[ERROR] No .fasta files found in {args.fasta_dir}")

    print(f"[INFO] Found {len(fasta_paths)} complexes to model.")

    for fa in fasta_paths:
        cid = Path(fa).stem
        out_dir = os.path.join(args.results_dir, cid)
        
        # Skip if already done (heuristic: check for done file or result existence)
        if os.path.exists(out_dir) and glob.glob(os.path.join(out_dir, "*.pdb")):
            print(f"[INFO] Skipping {cid}, results exist.")
            continue

        os.makedirs(out_dir, exist_ok=True)

        cmd = [
            "colabfold_batch",
            fa,
            out_dir,
            "--model-type", args.model_type,
            "--num-models", str(args.num_models),
            "--num-recycle", str(args.num_recycle),
            "--random-seed", str(args.seed),
        ]
        
        # Optional: Uncomment for AMBER relaxation (slower, better physical realism)
        # cmd += ["--amber"] 

        run_command(cmd, dry_run=args.dry_run)

    print(f"[SUCCESS] Batch processing complete. Results: {args.results_dir}")

if __name__ == "__main__":
    main()