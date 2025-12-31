#!/usr/bin/env python3
"""
AbDesign: Antibody Variant Library Generator
==================================================================

Companion script for the Part II paper:
"AI-powered integration of multi-source data for TAA discovery to accelerate 
ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization"

This script generates a library of CDR-focused antibody variants based on a 
benchmark antibody (Tezepelumab). It produces FASTA files formatted for 
structure prediction (Fv only) and docking (Ab:Ag complex).

Features:
---------
- CDR-focused mutagenesis preserving framework regions
- Generation of single-chain Fv (scFv) sequences with (GGGGS)3 linkers
- Preparation of multimer FASTA inputs for AlphaFold-Multimer
- Configurable mutation rates and random seeding for reproducibility

Requirements:
-------------
- python >= 3.8

Usage:
------
    python step0_generate_library.py --n_variants 100 --out_prefix tezepelumab_lib

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TADA
"""

import argparse
import random
from pathlib import Path

# ============================================================================
# Sequence Definitions (Tezepelumab & TSLP)
# ============================================================================

# Example: Tezepelumab Heavy Chain (VH) - Reference, replace it with your own sequences
HEAVY = ( 
    "QMQLVESGGGVVQPGRSLRLSCAASGFTFRTYGMHWVRQAPGKGLEWVAVIWYDGSNKHY"
    "ADSVKGRFTITRDNSKNTLNLQMNSLRAEDTAVYYCARAPQWELVHEAFDIWGQGTMVTV"
    "SS"
)

# Example: Tezepelumab Light Chain (VL) - Reference, replace it with your own sequences
LIGHT = ( 
    "SYVLTQPPSVSVAPGQTARITCGGNNLGSKSVHWYQQKPGQAPVLVVYDDSDRPSWIPER"
    "FSGSNSGNTATLTISRGEAGDEADYYCQVWDSSSDHVVFGGGTKLTVL"
)

VH_TEMPLATE = HEAVY[:122]
VL_TEMPLATE = LIGHT[:108]

# Example: TSLP Antigen Sequence (Uniprot Q969D9) - Reference, replace if needed
TSLP_SEQ = (
    "MFPFALLYVLSVSFRKIFILQLVGLVLTYDFTNCDFEKIKAAYLSTISKDLITYMSGTKS"
    "TEFNNTVSCSNRPHCLTEIQSLTFNPTAGCASLAKEMFAMKTKAALAIWCPGYSETQINA"
    "TQAMKKRRKRKVTTNKCLEQVSQLQGLWRRFNRPLLKQQ"
)

# ============================================================================
# CDR Configuration (Kabat Numbering)
# ============================================================================

# Light Chain CDRs (0-based indices converted from Kabat)
L1 = list(range(23 - 1, 33))      
L2 = list(range(49 - 1, 55))     
L3 = list(range(88 - 1, 98))   

# Heavy Chain CDRs (0-based indices)
H1 = list(range(26 - 1, 35))      
H2 = list(range(50 - 1, 66))      
H3 = list(range(98 - 1, 111))      

CDR_POS_H = sorted(set(H1 + H2 + H3))
CDR_POS_L = sorted(set(L1 + L2 + L3))

# Amino acid alphabet (excluding Cysteine to avoid unpaired disulfides)
AA_SET = "ACDEFGHIKLMNPQRSTVWY"  


# ============================================================================
# Mutagenesis Logic
# ============================================================================

def mutate_sequence(seq, positions, mut_fraction, rng=None):
    """
    Mutate a specific fraction of residues within defined positions.
    
    Args:
        seq (str): The amino acid sequence to mutate.
        positions (list): List of 0-based indices eligible for mutation.
        mut_fraction (float): Fraction of positions to mutate.
        rng (random.Random): Random number generator instance.
        
    Returns:
        tuple: (mutated_sequence, list_of_mutated_positions)
    """
    if rng is None:
        rng = random

    n_pos = len(positions)
    # Fixed count override based on paper methodology (example: 6 mutations)
    n_mut = 6 

    mut_positions = rng.sample(positions, n_mut)
    seq_list = list(seq)

    # print(f"[DEBUG] Mutating positions: {mut_positions}")

    for idx in mut_positions:
        orig = seq_list[idx]
        # Exclude Cysteine and the original amino acid
        allowed = [aa for aa in AA_SET if aa != orig and aa != "C"]
        new_aa = rng.choice(allowed)
        seq_list[idx] = new_aa

    return "".join(seq_list), mut_positions


def generate_library(n_variants: int, out_prefix: str, mut_fraction: float = 0.02, seed: int = 42):
    """
    Generate the variant library and write FASTA files.

    Args:
        n_variants (int): Number of variants to generate.
        out_prefix (str): Prefix for output filenames.
        mut_fraction (float): Target mutation rate.
        seed (int): Random seed for reproducibility.
    """
    rng = random.Random(seed)

    # Output paths
    fv_fasta = Path(f"{out_prefix}_fv.fasta")           # For structural checks
    complex_fasta = Path(f"{out_prefix}_complex.fasta") # For docking (AlphaFold)

    print(f"[INFO] Generating {n_variants} variants with seed {seed}...")

    with fv_fasta.open("w") as f_fv, complex_fasta.open("w") as f_cx:
        # Write reference (var0000)
        # Note: Logic implies 1-N generated, 0 is usually added manually or implied.
        
        for i in range(1, n_variants + 1):
            vh_mut, h_mutpos = mutate_sequence(VH_TEMPLATE, CDR_POS_H, mut_fraction, rng)
            vl_mut, l_mutpos = mutate_sequence(VL_TEMPLATE, CDR_POS_L, mut_fraction, rng)

            # Construct scFv: VH - (GGGGS)3 - VL
            linker = "GGGGS" * 3
            fv_seq = vh_mut + linker + vl_mut

            var_name = f"tezepelumab_var_{i:04d}"

            # 1. Fv FASTA (Single chain)
            f_fv.write(f">{var_name}_fv\n")
            f_fv.write(fv_seq + "\n")

            # 2. Complex FASTA (ColabFold multimer format: separated by colon)
            f_cx.write(f">{var_name}_complex\n")
            f_cx.write(f"{fv_seq}:{TSLP_SEQ}\n")

    print(f"[SUCCESS] Library generation complete.")
    print(f"   - Fv FASTA      : {fv_fasta}")
    print(f"   - Complex FASTA : {complex_fasta}")


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate CDR-mutated Tezepelumab variants as FASTA."
    )
    parser.add_argument("-n", "--n_variants", type=int, default=100, help="Number of variants (default: 100)")
    parser.add_argument("--prefix", type=str, default="tezepelumab_lib", help="Output prefix")
    parser.