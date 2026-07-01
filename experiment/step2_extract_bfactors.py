#!/usr/bin/env python3
"""
Step 2: Extract B-factors from binding-site residues.
======================================================
Usage:
    python step2_extract_bfactors.py

For each downloaded PDB structure:
    1. Identify the co-crystallized ligand
    2. Find binding-site residues within 5Å of the ligand
    3. Extract per-residue B-factors (average + CA only)
    4. Compute relative B-factor: B(residue) / B(pocket average)
    5. Save results to CSV

Output:
    output/bfactor_summary.csv   — per-residue B-factors for all structures
    output/bfactor_per_structure.csv — per-structure summary statistics
"""

import csv
import os
import sys

import numpy as np
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.dirname(__file__))
from config import PDB_DIR, OUTPUT_DIR, BINDING_SITE_DISTANCE

# Standard amino acid 3-letter codes
AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
}

# Common crystallization additives / ions to ignore
COMMON_LIGANDS = {
    'HOH', 'WAT', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'DMS',
    'ACT', 'ACE', 'NH2', 'NA', 'CL', 'MG', 'CA', 'ZN', 'FE',
    'K', 'BR', 'IOD', 'FMT', 'BME', 'DTT', 'TRS', 'HEPES',
}


def identify_ligand(pdb_path: str) -> str | None:
    """
    Identify the most likely co-crystallized drug ligand.
    Heuristic: largest non-standard, non-solvent, non-ion residue.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]

    ligand_sizes = {}  # resname -> number of heavy atoms
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in AA3 or resname in COMMON_LIGANDS:
                continue
            heavy_atoms = sum(1 for a in residue if a.element != 'H')
            if heavy_atoms > 5:  # Skip tiny fragments
                ligand_sizes[resname] = max(ligand_sizes.get(resname, 0), heavy_atoms)

    if not ligand_sizes:
        return None
    # Return the largest ligand
    return max(ligand_sizes, key=ligand_sizes.get)


def extract_bfactors(pdb_path: str, ligand_resname: str,
                     distance_cutoff: float = 5.0) -> list[dict]:
    """
    Extract per-residue B-factors for binding-site residues.

    Returns list of dicts with:
        chain, resseq, resname, avg_bfactor, ca_bfactor, n_atoms
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]

    # Collect ligand atom coordinates
    lig_coords = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip() == ligand_resname:
                for atom in residue:
                    if atom.element != 'H':
                        lig_coords.append(atom.get_coord())

    if not lig_coords:
        return []
    lig_coords = np.array(lig_coords)

    # Find nearby protein residues
    results = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname not in AA3:
                continue

            # Distance check
            min_dist = float('inf')
            for atom in residue:
                dists = np.linalg.norm(lig_coords - atom.get_coord(), axis=1)
                d = np.min(dists)
                if d < min_dist:
                    min_dist = d

            if min_dist >= distance_cutoff:
                continue

            # Extract B-factors
            bfactors = [a.get_bfactor() for a in residue]
            ca_bfactors = [a.get_bfactor() for a in residue if a.get_id() == "CA"]

            results.append({
                "chain": chain.get_id(),
                "resseq": residue.get_id()[1],
                "icode": residue.get_id()[2].strip(),
                "resname": resname,
                "avg_bfactor": round(np.mean(bfactors), 2),
                "ca_bfactor": round(ca_bfactors[0], 2) if ca_bfactors else None,
                "min_dist_to_ligand": round(min_dist, 2),
                "n_atoms": len(residue),
            })

    return results


def compute_relative_bfactors(residues: list[dict]) -> list[dict]:
    """
    Compute relative B-factor for each residue:
        relative_B = B(residue) / B(pocket_average)
    Analogous to BFIbs (Halip et al. 2021) but at residue level.
    """
    if not residues:
        return residues

    pocket_avg = np.mean([r["avg_bfactor"] for r in residues])
    pocket_ca_avg = np.mean([r["ca_bfactor"] for r in residues
                             if r["ca_bfactor"] is not None])

    for r in residues:
        r["relative_bfactor"] = round(r["avg_bfactor"] / pocket_avg, 3) if pocket_avg > 0 else 0
        r["relative_ca_bfactor"] = round(r["ca_bfactor"] / pocket_ca_avg, 3) \
            if r["ca_bfactor"] is not None and pocket_ca_avg > 0 else None
        r["pocket_avg_bfactor"] = round(pocket_avg, 2)
        r["pocket_ca_avg_bfactor"] = round(pocket_ca_avg, 2) if pocket_ca_avg else None

    return residues


def main():
    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
    if not pdb_files:
        print("ERROR: No PDB files found in data/pdb_files/")
        print("Run step1_download_structures.py --download first.")
        sys.exit(1)

    print(f"Processing {len(pdb_files)} structures...\n")

    all_residues = []
    structure_summary = []

    for pdb_file in pdb_files:
        pdb_id = pdb_file.replace(".pdb", "").upper()
        pdb_path = os.path.join(PDB_DIR, pdb_file)

        # Identify ligand
        ligand = identify_ligand(pdb_path)
        if not ligand:
            print(f"  [skip] {pdb_id}: no ligand found")
            continue

        # Extract binding-site B-factors
        residues = extract_bfactors(pdb_path, ligand, BINDING_SITE_DISTANCE)
        if not residues:
            print(f"  [skip] {pdb_id}: no binding-site residues within {BINDING_SITE_DISTANCE}Å")
            continue

        # Compute relative B-factors
        residues = compute_relative_bfactors(residues)

        # Add PDB ID and ligand to each residue record
        for r in residues:
            r["pdb_id"] = pdb_id
            r["ligand_resname"] = ligand

        all_residues.extend(residues)

        # Per-structure summary
        avg_b = np.mean([r["avg_bfactor"] for r in residues])
        std_b = np.std([r["avg_bfactor"] for r in residues])
        structure_summary.append({
            "pdb_id": pdb_id,
            "ligand": ligand,
            "n_binding_site_residues": len(residues),
            "pocket_avg_bfactor": round(avg_b, 2),
            "pocket_std_bfactor": round(std_b, 2),
            "max_bfactor_residue": max(residues, key=lambda r: r["avg_bfactor"])["resname"]
                                   + str(max(residues, key=lambda r: r["avg_bfactor"])["resseq"]),
        })
        print(f"  [ok] {pdb_id}: {len(residues)} residues, ligand={ligand}, "
              f"pocket_avg_B={avg_b:.1f}")

    # === Save per-residue CSV ===
    residue_csv = os.path.join(OUTPUT_DIR, "bfactor_summary.csv")
    fieldnames = ["pdb_id", "chain", "resseq", "resname", "avg_bfactor", "ca_bfactor",
                  "relative_bfactor", "relative_ca_bfactor", "pocket_avg_bfactor",
                  "min_dist_to_ligand", "ligand_resname"]
    with open(residue_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_residues)
    print(f"\nSaved per-residue data: {residue_csv} ({len(all_residues)} rows)")

    # === Save per-structure CSV ===
    struct_csv = os.path.join(OUTPUT_DIR, "bfactor_per_structure.csv")
    struct_fields = ["pdb_id", "ligand", "n_binding_site_residues",
                     "pocket_avg_bfactor", "pocket_std_bfactor", "max_bfactor_residue"]
    with open(struct_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=struct_fields)
        writer.writeheader()
        writer.writerows(structure_summary)
    print(f"Saved per-structure summary: {struct_csv} ({len(structure_summary)} rows)")


if __name__ == "__main__":
    main()
