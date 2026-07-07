#!/usr/bin/env python3
"""
Step 6: EnOpt — ML-based ensemble selection benchmark.
=======================================================
Usage:
    python step6_enopt.py --export      # Export EnOpt-compatible score matrix
    python step6_enopt.py --run         # Run EnOpt and parse ensemble selection
    python step6_enopt.py --all         # Export + run

This script bridges BFIbs-Ensemble with EnOpt (Bhatt et al., 2024):
    1. Export cross-docking affinity matrix in EnOpt CSV format
       (rows = ligands, columns = receptor conformations)
    2. Create known-ligands file (random subset of co-crystallized ligands)
    3. Run EnOpt to identify the most predictive sub-ensemble via ML
    4. Parse EnOpt output → selected PDB IDs for each ensemble size
    5. Save to output/enopt_ensemble.csv for step5 integration

Prerequisites:
    - Step 3 (cross-docking) must be completed
    - EnOpt cloned at experiment/enopt/
    - pip install pandas numpy scipy scikit-learn xgboost plotly
"""

import argparse
import csv
import os
import random
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import (
    OUTPUT_DIR, ENSEMBLE_SIZES,
)

ENOPT_DIR = os.path.join(os.path.dirname(__file__), "enopt")
ENOPT_SCRIPT = os.path.join(ENOPT_DIR, "ensemble_optimizer.py")


# =============================================================================
# Part 1: Export EnOpt-compatible score matrix
# =============================================================================

def export_enopt_matrix() -> str:
    """
    Convert cross_docking_results.csv to EnOpt format.
    
    EnOpt expects:
        - CSV with first column = ligand names, rest = conformation (receptor) columns
        - Values = docking scores (more negative = better for Vina)
    
    Returns:
        Path to the exported matrix CSV.
    """
    crossdock_file = os.path.join(OUTPUT_DIR, "cross_docking_results.csv")
    if not os.path.exists(crossdock_file):
        print("ERROR: cross_docking_results.csv not found. Run step3 first.")
        sys.exit(1)

    # Load cross-docking results
    rows = []
    with open(crossdock_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['affinity'] and row['affinity'] != 'None':
                rows.append({
                    'ligand': row['ligand_from'],
                    'receptor': row['receptor'],
                    'affinity': float(row['affinity']),
                })

    if not rows:
        print("ERROR: No valid docking results with affinity scores.")
        sys.exit(1)

    df = pd.DataFrame(rows)
    
    # Pivot: rows = ligands, columns = receptors, values = affinity
    matrix = df.pivot_table(
        index='ligand', columns='receptor', values='affinity', aggfunc='min'
    )
    
    # Fill self-docking diagonal (not computed in cross-docking) with best affinity
    # Rationale: native ligand should bind well to its own receptor
    best_affinity = matrix.min().min()  # most negative = best
    for pdb_id in matrix.index:
        if pdb_id in matrix.columns:
            if pd.isna(matrix.loc[pdb_id, pdb_id]):
                matrix.loc[pdb_id, pdb_id] = best_affinity
    
    # Fill any remaining NaN with 0 (poor binding proxy)
    matrix = matrix.fillna(0)
    
    # EnOpt expects first column = ligand names
    matrix = matrix.reset_index()
    matrix = matrix.rename(columns={'index': 'Ligand'})
    # Use original ligand name (keep PDB ID as ligand identifier)
    matrix['Ligand'] = matrix['Ligand'].apply(lambda x: f"ligand_{x}")
    
    # Save
    matrix_path = os.path.join(OUTPUT_DIR, "enopt_score_matrix.csv")
    matrix.to_csv(matrix_path, index=False)
    
    print(f"Exported EnOpt score matrix: {matrix_path}")
    print(f"  Shape: {matrix.shape[0]} ligands × {matrix.shape[1] - 1} conformations")
    return matrix_path


def create_known_ligands(matrix_path: str) -> str:
    """
    Create known-ligands file for EnOpt.
    
    Randomly selects ~50% of ligands as "known actives."
    All our ligands are co-crystallized CDK2 binders, so any split is valid.
    
    Returns:
        Path to the known ligands CSV.
    """
    matrix = pd.read_csv(matrix_path)
    all_ligands = matrix['Ligand'].tolist()
    
    # Random 50% split (seeded for reproducibility)
    random.seed(42)
    n_known = max(1, len(all_ligands) // 2)
    known_ligands = random.sample(all_ligands, n_known)
    
    known_path = os.path.join(OUTPUT_DIR, "enopt_known_ligands.csv")
    with open(known_path, 'w', newline='') as f:
        for lig in known_ligands:
            f.write(f"{lig}\n")
    
    print(f"Created known ligands file: {known_path}")
    print(f"  {len(known_ligands)} known / {len(all_ligands)} total ligands")
    return known_path


# =============================================================================
# Part 2: Run EnOpt
# =============================================================================

def run_enopt(matrix_path: str, known_path: str) -> str:
    """
    Run EnOpt as subprocess.
    
    Args:
        matrix_path: Path to EnOpt-compatible score matrix CSV.
        known_path: Path to known ligands CSV.
    
    Returns:
        Output file prefix used by EnOpt.
    """
    if not os.path.exists(ENOPT_SCRIPT):
        print(f"ERROR: EnOpt not found at {ENOPT_SCRIPT}")
        print("Clone it: git clone https://github.com/durrantlab/EnOpt.git experiment/enopt")
        sys.exit(1)
    
    out_prefix = os.path.join(OUTPUT_DIR, "enopt_output")
    max_size = max(ENSEMBLE_SIZES)  # Run with max ensemble size
    
    cmd = [
        sys.executable, ENOPT_SCRIPT,
        "-f", matrix_path,
        "-l", known_path,
        "--out_file", out_prefix,
        "--weighted_score",
        "--topn_confs", str(max_size),
    ]
    
    print(f"\nRunning EnOpt (topn_confs={max_size})...")
    print(f"  Command: {' '.join(cmd)}")
    
    result = subprocess.run(
        cmd,
        cwd=ENOPT_DIR,
        capture_output=True,
        text=True,
        timeout=300,  # 5 min timeout
    )
    
    if result.returncode != 0:
        print(f"EnOpt failed (exit code {result.returncode}):")
        print(result.stderr[-500:])
        sys.exit(1)
    
    print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
    print(f"EnOpt completed. Output prefix: {out_prefix}")
    return out_prefix


# =============================================================================
# Part 3: Parse EnOpt output → ensemble selection
# =============================================================================

def parse_enopt_ensemble(out_prefix: str) -> dict:
    """
    Parse EnOpt's conformations output to get the best ensemble.
    
    EnOpt outputs {prefix}_conformations.csv with:
        - "Conformation weight Model N" rows: feature importance per conformation
        - "Best subens. Model N" rows: True/False for top N conformations
    
    We aggregate across 3 CV models: sum the "Best subens." booleans,
    then rank conformations by total votes.
    
    Returns:
        {size: [pdb_ids]} — ensemble selection for each size in ENSEMBLE_SIZES.
    """
    conf_file = f"{out_prefix}_conformations.csv"
    if not os.path.exists(conf_file):
        print(f"ERROR: EnOpt output not found: {conf_file}")
        sys.exit(1)
    
    conf_df = pd.read_csv(conf_file, index_col=0)
    
    # Extract Best subens. rows (rows 1, 3, 5 → index 1, 3, 5)
    best_rows = conf_df.iloc[[1, 3, 5]]  # Best subens. Model 1, 2, 3
    
    # Sum True votes across models (True=1, False=0)
    vote_counts = best_rows.sum(axis=0)
    
    # Rank conformations by votes (descending), tie-break by weight sum
    weight_rows = conf_df.iloc[[0, 2, 4]]  # Conformation weight Model 1, 2, 3
    weight_sums = weight_rows.sum(axis=0)
    
    # Sort: first by votes (desc), then by weight sum (desc)
    ranking = pd.DataFrame({
        'conformation': vote_counts.index,
        'votes': vote_counts.values,
        'weight_sum': weight_sums.values,
    })
    ranking = ranking.sort_values(['votes', 'weight_sum'], ascending=[False, False])
    
    print("\nEnOpt conformation ranking (votes across 3 CV models):")
    for _, row in ranking.iterrows():
        bar = '█' * int(row['votes'])
        print(f"  {row['conformation']}: {int(row['votes'])} votes  {bar}")
    
    # Build ensemble selections for each size
    ensemble = {}
    all_confs = ranking['conformation'].tolist()
    
    for size in ENSEMBLE_SIZES:
        if size <= len(all_confs):
            ensemble[size] = [all_confs[:size]]
            print(f"  → Ensemble n={size}: {all_confs[:size]}")
    
    # Save to CSV for step5
    csv_path = os.path.join(OUTPUT_DIR, "enopt_ensemble.csv")
    rows = []
    for size, selections in ensemble.items():
        for i, sel in enumerate(selections):
            rows.append({
                'ensemble_size': size,
                'ensemble_id': i,
                'pdb_ids': ';'.join(sel),
                'strategy': 'enopt',
            })
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['strategy', 'ensemble_size', 'ensemble_id', 'pdb_ids'])
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"\nSaved EnOpt ensemble selections to: {csv_path}")
    return ensemble


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export', action='store_true', help='Export EnOpt score matrix')
    parser.add_argument('--run', action='store_true', help='Run EnOpt and parse ensemble')
    parser.add_argument('--all', action='store_true', help='Export + run')
    
    args = parser.parse_args()
    
    if not any([args.export, args.run, args.all]):
        parser.print_help()
        return
    
    if args.export or args.all:
        print("=" * 60)
        print("STEP 6a: Export EnOpt score matrix")
        print("=" * 60)
        matrix_path = export_enopt_matrix()
        known_path = create_known_ligands(matrix_path)
    
    if args.run or args.all:
        print("\n" + "=" * 60)
        print("STEP 6b: Run EnOpt")
        print("=" * 60)
        
        # Reload paths
        matrix_path = os.path.join(OUTPUT_DIR, "enopt_score_matrix.csv")
        known_path = os.path.join(OUTPUT_DIR, "enopt_known_ligands.csv")
        
        if not os.path.exists(matrix_path):
            print("Matrix not found. Run --export first.")
            return
        
        out_prefix = run_enopt(matrix_path, known_path)
        
        print("\n" + "=" * 60)
        print("STEP 6c: Parse EnOpt ensemble")
        print("=" * 60)
        parse_enopt_ensemble(out_prefix)


if __name__ == "__main__":
    main()
