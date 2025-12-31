#!/usr/bin/env python3
"""
AbDesign: AlphaFold2 Batch Runner
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script wraps ColabFold to execute batch structure prediction for the 
generated antibody library. It handles directory creation, query parsing, 
and execution of the AlphaFold2/AlphaFold-Multimer models.

Features:
---------
- Automated batch processing of FASTA files
- Integration with localColabFold
- Support for both monomer (Fv) and multimer (Ab:Ag) modes
- Automatic summary generation post-run

Requirements:
-------------
- colabfold
- localcolabfold environment

Usage:
------
    python step1.1_batch_af2_library.py input.fasta --outdir results/ --multimer

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import csv
import glob
import json
import os
from pathlib import Path

from colabfold.batch import get_queries, run

# ============================================================================
# Core Logic
# ============================================================================

def run_colabfold(fasta_path: str, outdir: str, model_type: str, use_multimer: bool):
    """
    Execute ColabFold batch run.
    """
    queries, is_complex = get_queries(fasta_path)
    result_dir = Path(outdir)
    result_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Loaded {len(queries)} queries from {fasta_path}")
    print(f"[INFO] is_complex (from FASTA): {is_complex}")
    print(f"[INFO] Running model_type      : {model_type}")

    run(
        queries=queries,
        result_dir=str(result_dir),
        use_amber=False,        # Disable AMBER relaxation for speed
        num_recycles=3,
        model_type=model_type,
        num_models=5,           # Generate 5 models per variant
        is_complex=use_multimer or is_complex,
    )
    print("[INFO] ColabFold run finished.")


def summarize_results(outdir: str, csv_path: str):
    """
    Generate a quick summary CSV of the Rank 1 models immediately after run.
    """
    outdir = Path(outdir)
    rows = []

    # Pattern for rank 001 JSON output from ColabFold
    pattern = str(outdir / "*_scores_rank_001_*.json")
    score_files = sorted(glob.glob(pattern))

    if not score_files:
        print(f"[WARN] No score json files found in {outdir}")
        return

    for jf in score_files:
        name = os.path.basename(jf)
        # Parse jobname: jobname_scores_rank_001_...
        jobname = name.split("_scores_rank_001_")[0]

        with open(jf, "r") as f:
            data = json.load(f)

        row = {
            "jobname": jobname,
            "mean_plddt": data.get("mean_plddt"),
            "ptm": data.get("ptm"),
            "iptm_plus_ptm": data.get("iptm+ptm"),
            "ranking_confidence": data.get("ranking_confidence"),
        }
        rows.append(row)

    # Write CSV
    fieldnames = ["jobname", "mean_plddt", "ptm", "iptm_plus_ptm", "ranking_confidence"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INFO] Wrote summary CSV: {csv_path} ({len(rows)} rows)")


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold on a library FASTA and summarize top-1 scores."
    )
    parser.add_argument("fasta", help="Input FASTA file (Fv or complex)")
    parser.add_argument(
        "-o", "--outdir", 
        default="af2_results",
        help="Output directory for ColabFold results (default: af2_results)"
    )
    parser.add_argument(
        "--model_type", 
        default="alphafold2_ptm",
        help="Model type (alphafold2_ptm or alphafold2_multimer_v3, etc.)"
    )
    parser.add_argument(
        "--multimer", 
        action="store_true",
        help="Force multimer mode (necessary for Ab:Ag complexes)"
    )
    args = parser.parse_args()

    run_colabfold(args.fasta, args.outdir, args.model_type, args.multimer)

    # Post-run summary
    csv_path = Path(args.outdir) / "af2_summary_rank001.csv"
    summarize_results(args.outdir, str(csv_path))


if __name__ == "__main__":
    main()