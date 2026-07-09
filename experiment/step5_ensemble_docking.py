#!/usr/bin/env python3
"""
Step 5: B-factor guided ensemble docking.
==========================================
Usage:
    python step5_ensemble_docking.py --select     # Select ensembles by different strategies
    python step5_ensemble_docking.py --dock       # Run ensemble docking
    python step5_ensemble_docking.py --compare    # Compare strategies
    python step5_ensemble_docking.py --all        # Run everything

This script:
    1. Calculate BFIbs (B-factor Index for binding site) for each structure
    2. Select ensemble members using different strategies:
       - B-factor guided: BFIbs closest to 1 (balanced flexibility)
       - Random: baseline comparison
       - Lowest B-factor: most rigid binding sites
    3. Run ensemble docking: for each ligand, dock against all ensemble members,
       take the best pose (lowest RMSD to crystal)
    4. Compare success rates across strategies

Output:
    output/bfibs_scores.csv              — BFIbs for each structure
    output/ensemble_selections.csv       — Which structures in each ensemble
    output/ensemble_docking_results.csv  — Per-ligand ensemble docking results
    output/ensemble_comparison.csv       — Strategy comparison summary
    figures/ensemble_comparison.png      — Bar chart comparing strategies
"""

import argparse
import csv
import os
import random
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.dirname(__file__))
from config import (
    PDB_DIR, OUTPUT_DIR, RESULTS_DIR, FIGURES_DIR,
    BINDING_SITE_DISTANCE, ENSEMBLE_SIZES, ENSEMBLE_STRATEGIES,
    NUM_RANDOM_ENSEMBLES, RMSD_SUCCESS_THRESHOLD
)

# Standard amino acid 3-letter codes
AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
}


# =============================================================================
# Part 1: Calculate BFIbs for each structure
# =============================================================================

def calculate_bfibs(pdb_id: str) -> dict:
    """
    Calculate BFIbs = median(B_pocket) / median(B_protein)
    
    BFIbs < 1: binding site is more rigid than average
    BFIbs ≈ 1: binding site flexibility matches protein average
    BFIbs > 1: binding site is more flexible than average
    """
    # Try both lowercase and uppercase filenames
    for name in (pdb_id.lower(), pdb_id.upper()):
        pdb_file = os.path.join(PDB_DIR, f"{name}.pdb")
        if os.path.exists(pdb_file):
            break
    else:
        return None
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]
    
    # Find ligand (first non-protein, non-water HETATM)
    ligand_atoms = []
    protein_atoms = []
    
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in ('HOH', 'WAT', 'DOD'):
                continue
            if resname in AA3:
                for atom in residue:
                    if atom.get_name() != 'H':
                        protein_atoms.append(atom)
            elif len(resname.strip()) > 1:
                for atom in residue:
                    if atom.get_name() != 'H':
                        ligand_atoms.append(atom)
    
    if not ligand_atoms or not protein_atoms:
        return None
    
    # Get ligand heavy atom coordinates
    ligand_coords = np.array([a.get_coord() for a in ligand_atoms])
    
    # Find binding site residues
    pocket_bfactors = []
    protein_bfactors = []
    
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname not in AA3:
                continue
            
            # Calculate distance to ligand
            min_dist = float('inf')
            for atom in residue:
                if atom.get_name() == 'H':
                    continue
                coord = atom.get_coord()
                for lig_coord in ligand_coords:
                    d = np.linalg.norm(coord - lig_coord)
                    min_dist = min(min_dist, d)
            
            # Get residue B-factors
            res_bfactors = [a.get_bfactor() for a in residue 
                          if a.get_name() != 'H']
            
            if not res_bfactors:
                continue
            
            protein_bfactors.extend(res_bfactors)
            
            if min_dist <= BINDING_SITE_DISTANCE:
                pocket_bfactors.extend(res_bfactors)
    
    if not pocket_bfactors or not protein_bfactors:
        return None
    
    # Calculate BFIbs
    median_pocket = np.median(pocket_bfactors)
    median_protein = np.median(protein_bfactors)
    
    if median_protein == 0:
        return None
    
    bfibs = median_pocket / median_protein
    
    return {
        'pdb_id': pdb_id,
        'bfibs': bfibs,
        'median_pocket_bfactor': median_pocket,
        'median_protein_bfactor': median_protein,
        'mean_pocket_bfactor': np.mean(pocket_bfactors),
        'num_pocket_residues': len(pocket_bfactors),
    }


def compute_all_bfibs() -> list:
    """Calculate BFIbs for all downloaded structures."""
    pdb_files = [f for f in os.listdir(PDB_DIR) if f.endswith('.pdb')]
    # Normalize to uppercase for consistency with cross-docking results
    pdb_ids = [f.replace('.pdb', '').upper() for f in pdb_files]
    
    results = []
    for pdb_id in sorted(pdb_ids):
        result = calculate_bfibs(pdb_id)
        if result:
            results.append(result)
            print(f"  {pdb_id}: BFIbs = {result['bfibs']:.3f}")
    
    # Save to CSV
    if results:
        csv_file = os.path.join(OUTPUT_DIR, "bfibs_scores.csv")
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSaved BFIbs scores to: {csv_file}")
    
    return results


# =============================================================================
# Part 2: Ensemble selection strategies
# =============================================================================

def select_ensemble_bfactor_guided(bfibs_data: list, size: int) -> list:
    """
    Select structures with BFIbs closest to 1.0
    Rationale: BFIbs ≈ 1 means binding site flexibility matches protein average,
    which should give the most "representative" conformations for docking.
    """
    sorted_data = sorted(bfibs_data, key=lambda x: abs(x['bfibs'] - 1.0))
    return [d['pdb_id'] for d in sorted_data[:size]]


def select_ensemble_lowest_bfactor(bfibs_data: list, size: int) -> list:
    """
    Select structures with lowest BFIbs (most rigid binding sites).
    Rationale: rigid binding sites should be easier for rigid docking.
    """
    sorted_data = sorted(bfibs_data, key=lambda x: x['bfibs'])
    return [d['pdb_id'] for d in sorted_data[:size]]


def select_ensemble_random(bfibs_data: list, size: int) -> list:
    """Random selection — baseline for comparison."""
    return random.sample([d['pdb_id'] for d in bfibs_data], 
                         min(size, len(bfibs_data)))


def select_all_ensembles(bfibs_data: list) -> dict:
    """
    Generate ensembles for all strategies and sizes.
    Returns: {strategy: {size: [pdb_ids]}}
    """
    ensembles = {}
    
    for strategy in ENSEMBLE_STRATEGIES:
        ensembles[strategy] = {}
        
        for size in ENSEMBLE_SIZES:
            if strategy == "enopt":
                # EnOpt ensemble is pre-computed by step6; read from CSV
                enopt_file = os.path.join(OUTPUT_DIR, "enopt_ensemble.csv")
                if not os.path.exists(enopt_file):
                    print(f"  [skip] enopt: {enopt_file} not found. Run step6 first.")
                    continue
                ensembles[strategy][size] = _load_enopt_ensemble(enopt_file, size)
                if not ensembles[strategy][size]:
                    del ensembles[strategy][size]
                continue
            
            if size > len(bfibs_data):
                continue
            
            if strategy == "bfactor_guided":
                ensembles[strategy][size] = [
                    select_ensemble_bfactor_guided(bfibs_data, size)
                ]
            elif strategy == "lowest_bfactor":
                ensembles[strategy][size] = [
                    select_ensemble_lowest_bfactor(bfibs_data, size)
                ]
            elif strategy == "random":
                # Generate multiple random ensembles for statistical comparison
                ensembles[strategy][size] = []
                for _ in range(NUM_RANDOM_ENSEMBLES):
                    ens = select_ensemble_random(bfibs_data, size)
                    ensembles[strategy][size].append(ens)
    
    return ensembles


def _load_enopt_ensemble(enopt_file: str, size: int) -> list:
    """Load EnOpt-selected ensemble for a given size from CSV."""
    selections = []
    with open(enopt_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row['ensemble_size']) == size:
                pdb_ids = row['pdb_ids'].split(';')
                selections.append(pdb_ids)
    return selections


def save_ensemble_selections(ensembles: dict, bfibs_data: list):
    """Save ensemble selections to CSV."""
    csv_file = os.path.join(OUTPUT_DIR, "ensemble_selections.csv")
    
    rows = []
    for strategy, sizes in ensembles.items():
        for size, ensemble_list in sizes.items():
            for i, ens in enumerate(ensemble_list):
                row = {
                    'strategy': strategy,
                    'size': size,
                    'ensemble_id': i,
                    'pdb_ids': ';'.join(ens),
                    'bfibs_values': ';'.join([
                        f"{next((d['bfibs'] for d in bfibs_data if d['pdb_id'] == pid), 0):.3f}"
                        for pid in ens
                    ])
                }
                rows.append(row)
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Saved ensemble selections to: {csv_file}")


# =============================================================================
# Part 3: Ensemble docking
# =============================================================================

def run_ensemble_docking(ensembles: dict):
    """
    For each ensemble, run docking and evaluate.
    
    Ensemble docking logic:
    - For each ligand (from structure j), dock against ALL receptors in ensemble
    - Take the BEST result (lowest RMSD to crystal pose of ligand j)
    - Success = best RMSD < threshold
    """
    # Load cross-docking results from step3
    crossdock_file = os.path.join(OUTPUT_DIR, "cross_docking_results.csv")
    if not os.path.exists(crossdock_file):
        print("ERROR: Run step3 cross-docking first!")
        return None
    
    # Load results into lookup: {(receptor, ligand): rmsd}
    results_lookup = {}
    with open(crossdock_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row.get('rmsd'):
                continue
            key = (row['receptor'], row['ligand_from'])
            results_lookup[key] = float(row['rmsd'])
    
    # Evaluate each ensemble
    all_results = []
    
    for strategy, sizes in ensembles.items():
        for size, ensemble_list in sizes.items():
            for ens_id, ens_pdb_ids in enumerate(ensemble_list):
                # For each ligand, find best RMSD across ensemble
                ligand_success = []
                ligand_best_rmsds = []
                
                # Get all ligands that are NOT in the ensemble
                all_ligands = set()
                for key in results_lookup:
                    all_ligands.add(key[1])
                
                for ligand in all_ligands:
                    best_rmsd = float('inf')
                    
                    for receptor in ens_pdb_ids:
                        key = (receptor, ligand)
                        if key in results_lookup:
                            rmsd = results_lookup[key]
                            best_rmsd = min(best_rmsd, rmsd)
                    
                    if best_rmsd < float('inf'):
                        success = 1 if best_rmsd < RMSD_SUCCESS_THRESHOLD else 0
                        ligand_success.append(success)
                        ligand_best_rmsds.append(best_rmsd)
                
                if ligand_success:
                    success_rate = np.mean(ligand_success)
                    mean_rmsd = np.mean(ligand_best_rmsds)
                    
                    all_results.append({
                        'strategy': strategy,
                        'ensemble_size': size,
                        'ensemble_id': ens_id,
                        'pdb_ids': ';'.join(ens_pdb_ids),
                        'success_rate': success_rate,
                        'mean_rmsd': mean_rmsd,
                        'median_rmsd': np.median(ligand_best_rmsds),
                        'num_ligands_tested': len(ligand_success),
                    })
    
    # Save results
    if all_results:
        csv_file = os.path.join(OUTPUT_DIR, "ensemble_docking_results.csv")
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
        print(f"Saved ensemble docking results to: {csv_file}")
    
    return all_results


# =============================================================================
# Part 4: Compare strategies
# =============================================================================

def compare_strategies(results: list):
    """Compare ensemble docking strategies."""
    if not results:
        print("No results to compare!")
        return
    
    # Group by strategy and size
    summary = defaultdict(list)
    for r in results:
        key = f"{r['strategy']}_n{r['ensemble_size']}"
        summary[key].append(r['success_rate'])
    
    print("\n" + "=" * 70)
    print("ENSEMBLE DOCKING COMPARISON")
    print("=" * 70)
    print(f"\n{'Strategy':<20} {'Size':<6} {'Success Rate':<15} {'Std':<8} {'N'}")
    print("-" * 70)
    
    comparison_rows = []
    for key in sorted(summary.keys()):
        rates = summary[key]
        strategy, size_str = key.rsplit('_n', 1)
        size = int(size_str)
        mean_rate = np.mean(rates)
        std_rate = np.std(rates) if len(rates) > 1 else 0
        
        print(f"{strategy:<20} {size:<6} {mean_rate:<15.3f} {std_rate:<8.3f} {len(rates)}")
        
        comparison_rows.append({
            'strategy': strategy,
            'ensemble_size': size,
            'mean_success_rate': mean_rate,
            'std_success_rate': std_rate,
            'num_ensembles': len(rates),
        })
    
    # Save comparison
    csv_file = os.path.join(OUTPUT_DIR, "ensemble_comparison.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=comparison_rows[0].keys())
        writer.writeheader()
        writer.writerows(comparison_rows)
    
    print(f"\nSaved comparison to: {csv_file}")
    
    # Generate plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        strategies = sorted(set(r['strategy'] for r in comparison_rows))
        colors = {'bfactor_guided': '#2ecc71', 'random': '#95a5a6', 
                  'lowest_bfactor': '#3498db', 'enopt': '#e74c3c'}
        
        x = np.arange(len(ENSEMBLE_SIZES))
        width = 0.25
        
        for i, strategy in enumerate(strategies):
            means = []
            stds = []
            for size in ENSEMBLE_SIZES:
                row = next((r for r in comparison_rows 
                           if r['strategy'] == strategy and r['ensemble_size'] == size), None)
                if row:
                    means.append(row['mean_success_rate'])
                    stds.append(row['std_success_rate'])
                else:
                    means.append(0)
                    stds.append(0)
            
            color = colors.get(strategy, '#e74c3c')
            ax.bar(x + i*width, means, width, label=strategy, 
                   yerr=stds, color=color, alpha=0.8, capsize=3)
        
        ax.set_xlabel('Ensemble Size')
        ax.set_ylabel('Success Rate (RMSD < 2Å)')
        ax.set_title('Ensemble Docking: B-factor Guided vs Random vs Lowest B-factor')
        ax.set_xticks(x + width)
        ax.set_xticklabels([str(s) for s in ENSEMBLE_SIZES])
        ax.legend()
        ax.set_ylim(0, 1)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "ensemble_comparison.png"), dpi=150)
        print(f"Saved figure to: {os.path.join(FIGURES_DIR, 'ensemble_comparison.png')}")
        
    except ImportError:
        print("matplotlib not available, skipping plot generation")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--select', action='store_true', help='Select ensembles')
    parser.add_argument('--dock', action='store_true', help='Run ensemble docking')
    parser.add_argument('--compare', action='store_true', help='Compare strategies')
    parser.add_argument('--all', action='store_true', help='Run everything')
    
    args = parser.parse_args()
    
    if not any([args.select, args.dock, args.compare, args.all]):
        parser.print_help()
        return
    
    # Step 1: Calculate BFIbs
    if args.select or args.all:
        print("=" * 60)
        print("Calculating BFIbs for all structures...")
        print("=" * 60)
        bfibs_data = compute_all_bfibs()
        
        if not bfibs_data:
            print("ERROR: No BFIbs data computed. Check PDB files.")
            return
        
        print("\nSelecting ensembles...")
        ensembles = select_all_ensembles(bfibs_data)
        save_ensemble_selections(ensembles, bfibs_data)
    
    # Step 2: Run ensemble docking
    if args.dock or args.all:
        print("\n" + "=" * 60)
        print("Running ensemble docking...")
        print("=" * 60)
        
        # Reload ensembles
        bfibs_data = compute_all_bfibs()
        ensembles = select_all_ensembles(bfibs_data)
        results = run_ensemble_docking(ensembles)
    
    # Step 3: Compare
    if args.compare or args.all:
        print("\n" + "=" * 60)
        print("Comparing strategies...")
        print("=" * 60)
        
        # Load results
        csv_file = os.path.join(OUTPUT_DIR, "ensemble_docking_results.csv")
        if os.path.exists(csv_file):
            results = []
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    row['ensemble_size'] = int(row['ensemble_size'])
                    row['success_rate'] = float(row['success_rate'])
                    row['mean_rmsd'] = float(row['mean_rmsd'])
                    row['median_rmsd'] = float(row['median_rmsd'])
                    row['num_ligands_tested'] = int(row['num_ligands_tested'])
                    results.append(row)
            compare_strategies(results)
        else:
            print("No ensemble docking results found. Run --dock first.")


if __name__ == "__main__":
    main()
