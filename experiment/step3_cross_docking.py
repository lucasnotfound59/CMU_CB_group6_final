#!/usr/bin/env python3
"""
Step 3: Prepare and run cross-docking with AutoDock Vina.
==========================================================
Usage:
    python step3_cross_docking.py --prepare    # Generate PDBQT files
    python step3_cross_docking.py --run        # Run Vina cross-docking
    python step3_cross_docking.py --all        # Both

Pipeline:
    For each pair of structures (receptor_i, ligand_j) where i ≠ j:
        1. Prepare receptor PDBQT (from structure i, remove ligand)
        2. Prepare ligand PDBQT (from structure j)
        3. Run Vina docking (ligand j → receptor i)
        4. Record best pose RMSD vs crystal pose of ligand j in structure j

Prerequisites:
    - AutoDock Vina installed and accessible (see config.VINA_BINARY)
    - Open Babel (obabel) or MGLTools (prepare_receptor4.py) for PDBQT conversion
    - step1 and step2 completed
"""

import argparse
import csv
import os
import subprocess
import sys

import numpy as np
from Bio.PDB import PDBParser, PDBIO

sys.path.insert(0, os.path.dirname(__file__))
from config import (PDB_DIR, VINA_DIR, RESULTS_DIR, OUTPUT_DIR,
                    VINA_BINARY, VINA_BOX_SIZE, VINA_EXHAUSTIVENESS,
                    VINA_NUM_MODES, BINDING_SITE_DISTANCE)

AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
}
COMMON_LIGANDS = {
    'HOH', 'WAT', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'DMS',
    'ACT', 'ACE', 'NA', 'CL', 'MG', 'CA', 'ZN', 'K', 'BR',
}


# === Helper: prepare receptor (protein only) ===

def prepare_receptor_pdb(pdb_path: str, ligand_resname: str, output_path: str):
    """
    Strip all non-protein atoms from PDB, save clean receptor PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_path)

    io = PDBIO()
    io.set_structure(structure)

    class ProteinOnly:
        def accept_residue(self, residue):
            return residue.get_resname().strip() in AA3

    io.save(output_path, ProteinOnly())


def prepare_ligand_pdb(pdb_path: str, ligand_resname: str, output_path: str):
    """
    Extract only the ligand from a PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ligand", pdb_path)
    model = structure[0]

    serial = 0
    with open(output_path, "w") as f:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() != ligand_resname:
                    continue
                for atom in residue:
                    serial += 1
                    # Infer element from atom name if not set
                    elem = atom.element if atom.element and atom.element != ' ' else atom.get_name()[0]
                    line = (
                        f"HETATM{serial:>5d} {atom.get_name():<4s}"
                        f"{residue.get_resname():>3s} "
                        f"{chain.get_id()}{residue.get_id()[1]:>4d}    "
                        f"{atom.get_coord()[0]:>8.3f}"
                        f"{atom.get_coord()[1]:>8.3f}"
                        f"{atom.get_coord()[2]:>8.3f}"
                        f"  1.00  0.00          {elem:>2s}\n"
                    )
                    f.write(line)
                f.write("TER\n")
                return  # Only first matching residue


def pdb_to_pdbqt(pdb_path: str, pdbqt_path: str, is_receptor: bool = True):
    """
    Convert PDB → PDBQT using Open Babel or MGLTools.
    Tries obabel first, then prepare_receptor4/prepare_ligand4.
    """
    # Try Open Babel
    try:
        if is_receptor:
            cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "-h"]
        else:
            cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-h", "--partialcharge", "gasteiger"]
        subprocess.run(cmd, check=True, capture_output=True)
        if os.path.exists(pdbqt_path):
            return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass

    # Try MGLTools
    try:
        if is_receptor:
            script = "prepare_receptor4.py"
            cmd = ["python", script, "-r", pdb_path, "-o", pdbqt_path]
        else:
            script = "prepare_ligand4.py"
            cmd = ["python", script, "-l", pdb_path, "-o", pdbqt_path]
        subprocess.run(cmd, check=True, capture_output=True)
        if os.path.exists(pdbqt_path):
            return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass

    print(f"    [error] Cannot convert {pdb_path} → PDBQT")
    print(f"            Install Open Babel (obabel) or MGLTools")
    return False


# === Helper: compute binding site center ===

def get_binding_site_center(pdb_path: str, ligand_resname: str) -> tuple:
    """Get the geometric center of the ligand (for Vina box center)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]

    coords = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip() == ligand_resname:
                for atom in residue:
                    if atom.element != 'H':
                        coords.append(atom.get_coord())

    if not coords:
        return (0, 0, 0)
    center = np.mean(coords, axis=0)
    return tuple(center)


# === Helper: identify ligand ===

def identify_ligand(pdb_path: str) -> str | None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = structure[0]
    ligand_sizes = {}
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in AA3 or resname in COMMON_LIGANDS:
                continue
            heavy = sum(1 for a in residue if a.element != 'H')
            if heavy > 5:
                ligand_sizes[resname] = max(ligand_sizes.get(resname, 0), heavy)
    return max(ligand_sizes, key=ligand_sizes.get) if ligand_sizes else None


# === Vina docking ===

def run_vina(receptor_pdbqt: str, ligand_pdbqt: str, output_pdbqt: str,
             center: tuple, box_size: tuple) -> dict | None:
    """
    Run AutoDock Vina and return result dict.
    """
    cmd = [
        VINA_BINARY,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--exhaustiveness", str(VINA_EXHAUSTIVENESS),
        "--num_modes", str(VINA_NUM_MODES),
        "--out", output_pdbqt,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = result.stdout + result.stderr

        # Parse Vina output for affinity
        best_affinity = None
        for line in output.split("\n"):
            if line.strip().startswith("1"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        best_affinity = float(parts[1])
                    except ValueError:
                        pass

        return {"affinity": best_affinity, "output": output, "success": result.returncode == 0}
    except subprocess.TimeoutExpired:
        return {"affinity": None, "output": "TIMEOUT", "success": False}
    except FileNotFoundError:
        print(f"    [error] Vina not found at: {VINA_BINARY}")
        print(f"            Install AutoDock Vina: https://vina.scripps.edu/")
        return None


# === Compute RMSD between two PDB files ===

def compute_rmsd(ref_path: str, docked_path: str, ref_resname: str) -> float | None:
    """
    Compute heavy-atom RMSD between reference ligand and best docked pose.
    """
    parser = PDBParser(QUIET=True)

    try:
        ref_struct = parser.get_structure("ref", ref_path)
        dock_struct = parser.get_structure("dock", docked_path)
    except Exception:
        return None

    # Reference ligand coords
    ref_coords = []
    for model in ref_struct:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ref_resname:
                    for atom in residue:
                        if atom.element != 'H':
                            ref_coords.append(atom.get_coord())
                    break
            break

    # Docked ligand coords (first model, first chain, first residue)
    dock_coords = []
    for model in dock_struct:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname not in AA3 and resname not in COMMON_LIGANDS:
                    for atom in residue:
                        if atom.element != 'H':
                            dock_coords.append(atom.get_coord())
                    break
            break

    if not ref_coords or not dock_coords:
        return None

    ref_arr = np.array(ref_coords)
    dock_arr = np.array(dock_coords)

    # Must have same number of atoms
    n = min(len(ref_arr), len(dock_arr))
    if n == 0:
        return None

    diff = ref_arr[:n] - dock_arr[:n]
    rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return round(rmsd, 3)


# === Main ===

def prepare_all():
    """Prepare PDBQT files for all structures."""
    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
    print(f"Preparing PDBQT files for {len(pdb_files)} structures...\n")

    for pdb_file in pdb_files:
        pdb_id = pdb_file.replace(".pdb", "").upper()
        pdb_path = os.path.join(PDB_DIR, pdb_file)
        ligand = identify_ligand(pdb_path)
        if not ligand:
            print(f"  [skip] {pdb_id}: no ligand")
            continue

        # Receptor
        rec_pdb = os.path.join(VINA_DIR, f"{pdb_id.lower()}_receptor.pdb")
        rec_pdbqt = os.path.join(VINA_DIR, f"{pdb_id.lower()}_receptor.pdbqt")
        prepare_receptor_pdb(pdb_path, ligand, rec_pdb)
        ok = pdb_to_pdbqt(rec_pdb, rec_pdbqt, is_receptor=True)
        if ok:
            print(f"  [ok] {pdb_id} receptor → PDBQT")

        # Ligand
        lig_pdb = os.path.join(VINA_DIR, f"{pdb_id.lower()}_ligand.pdb")
        lig_pdbqt = os.path.join(VINA_DIR, f"{pdb_id.lower()}_ligand.pdbqt")
        prepare_ligand_pdb(pdb_path, ligand, lig_pdb)
        ok = pdb_to_pdbqt(lig_pdb, lig_pdbqt, is_receptor=False)
        if ok:
            print(f"  [ok] {pdb_id} ligand → PDBQT")


def run_cross_docking():
    """Run all cross-docking pairs."""
    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
    structures = []
    for pf in pdb_files:
        pdb_id = pf.replace(".pdb", "").upper()
        pdb_path = os.path.join(PDB_DIR, pf)
        lig = identify_ligand(pdb_path)
        if lig:
            structures.append({"pdb_id": pdb_id, "ligand": lig, "path": pdb_path})

    print(f"Cross-docking {len(structures)} structures "
          f"({len(structures) * (len(structures) - 1)} pairs)...\n")

    results = []
    total = len(structures) * (len(structures) - 1)
    count = 0

    for rec in structures:
        rec_pdbqt = os.path.join(VINA_DIR, f"{rec['pdb_id'].lower()}_receptor.pdbqt")
        if not os.path.exists(rec_pdbqt):
            continue

        for lig_struct in structures:
            if lig_struct["pdb_id"] == rec["pdb_id"]:
                continue  # Skip self-docking
            count += 1

            lig_pdbqt = os.path.join(VINA_DIR, f"{lig_struct['pdb_id'].lower()}_ligand.pdbqt")
            if not os.path.exists(lig_pdbqt):
                continue

            # Vina output
            out_pdbqt = os.path.join(RESULTS_DIR,
                                     f"{rec['pdb_id'].lower()}_{lig_struct['pdb_id'].lower()}.pdbqt")

            # Box center = ligand center in the REFERENCE structure
            center = get_binding_site_center(lig_struct["path"], lig_struct["ligand"])

            # Run Vina
            print(f"  [{count}/{total}] {lig_struct['pdb_id']} → {rec['pdb_id']} ...", end=" ")
            vina_result = run_vina(rec_pdbqt, lig_pdbqt, out_pdbqt,
                                   center, VINA_BOX_SIZE)

            if vina_result is None or not vina_result["success"]:
                print("FAILED")
                results.append({
                    "receptor": rec["pdb_id"],
                    "ligand_from": lig_struct["pdb_id"],
                    "rmsd": None,
                    "affinity": None,
                    "status": "vina_failed",
                })
                continue

            # Compute RMSD vs crystal pose
            rmsd = compute_rmsd(lig_struct["path"], out_pdbqt, lig_struct["ligand"])
            status = "success" if rmsd is not None and rmsd < 2.0 else "high_rmsd"
            print(f"RMSD={rmsd}, affinity={vina_result['affinity']}")

            results.append({
                "receptor": rec["pdb_id"],
                "ligand_from": lig_struct["pdb_id"],
                "rmsd": rmsd,
                "affinity": vina_result["affinity"],
                "status": status,
            })

    # Save results
    csv_path = os.path.join(OUTPUT_DIR, "cross_docking_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["receptor", "ligand_from", "rmsd",
                                               "affinity", "status"])
        writer.writeheader()
        writer.writerows(results)

    success = sum(1 for r in results if r["status"] == "success")
    print(f"\nResults: {csv_path}")
    print(f"Success (RMSD < 2Å): {success}/{len(results)} "
          f"({100*success/len(results):.1f}%)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true", help="Prepare PDBQT files")
    parser.add_argument("--run", action="store_true", help="Run cross-docking")
    parser.add_argument("--all", action="store_true", help="Prepare + run")
    args = parser.parse_args()

    if args.all:
        prepare_all()
        run_cross_docking()
    elif args.prepare:
        prepare_all()
    elif args.run:
        run_cross_docking()
    else:
        parser.print_help()
