#!/usr/bin/env python3
"""
Positive control: pocket Cα RMSD vs cross-docking RMSD.
=========================================================
Computes pairwise binding-site structural similarity and tests
whether it predicts cross-docking accuracy better than B-factor.

Usage:
    BFIBS_TARGET_NAME=CDK1 python step4b_pocket_rmsd.py
"""

import csv, os, sys
import numpy as np
from Bio.PDB import PDBParser
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from config import PDB_DIR, OUTPUT_DIR, FIGURES_DIR, BINDING_SITE_DISTANCE

AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
       'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
COMMON = {'HOH','WAT','SO4','PO4','GOL','EDO','PEG','DMS',
          'ACT','ACE','NA','CL','MG','CA','ZN','K','BR'}

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def identify_ligand(pdb_path):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", pdb_path)
    sizes = {}
    for chain in s[0]:
        for res in chain:
            rn = res.get_resname().strip()
            if rn in AA3 or rn in COMMON: continue
            heavy = sum(1 for a in res if a.element != 'H')
            if heavy > 5:
                sizes[rn] = max(sizes.get(rn, 0), heavy)
    return max(sizes, key=sizes.get) if sizes else None


def get_binding_site_ca(pdb_path, ligand_resname):
    """Return {resseq: Cα_coord} for binding site residues."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", pdb_path)
    model = s[0]

    # Ligand coords
    lig_coords = []
    for chain in model:
        for res in chain:
            if res.get_resname().strip() == ligand_resname:
                for atom in res:
                    lig_coords.append(atom.get_coord())
    if not lig_coords:
        return {}
    lig_arr = np.array(lig_coords)

    # Binding site Cα atoms
    ca_dict = {}
    for chain in model:
        for res in chain:
            if res.get_resname().strip() not in AA3:
                continue
            ca = res.child_dict.get("CA")
            if ca is None:
                continue
            # Check if any atom within cutoff
            res_coords = np.array([a.get_coord() for a in res])
            dists = np.min(np.linalg.norm(
                lig_arr[:, None, :] - res_coords[None, :, :], axis=2))
            if dists < BINDING_SITE_DISTANCE:
                ca_dict[res.get_id()[1]] = ca.get_coord()
    return ca_dict


def kabsch_rmsd(coords_a, coords_b):
    """Kabsch-aligned RMSD between two coordinate sets (same residue order)."""
    a, b = np.array(coords_a), np.array(coords_b)
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    ac, bc = a - ca, b - cb
    H = bc.T @ ac
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1; R = Vt.T @ U.T
    return np.sqrt(np.mean(np.sum((ac - bc @ R.T) ** 2, axis=1)))


def compute_pairwise_pocket_rmsd():
    """Compute pocket Cα RMSD for all pairs of structures."""
    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
    structures = {}
    for pf in pdb_files:
        pdb_id = pf.replace(".pdb", "").upper()
        pdb_path = os.path.join(PDB_DIR, pf)
        lig = identify_ligand(pdb_path)
        if lig:
            structures[pdb_id] = get_binding_site_ca(pdb_path, lig)

    ids = sorted(structures.keys())
    pairs = []
    for rec_id in ids:
        for lig_id in ids:
            if rec_id == lig_id:
                continue
            ca_rec = structures[rec_id]
            ca_lig = structures[lig_id]
            # Align on common binding site residues
            common = sorted(set(ca_rec.keys()) & set(ca_lig.keys()))
            if len(common) < 5:
                continue
            rmsd = kabsch_rmsd(
                [ca_rec[k] for k in common],
                [ca_lig[k] for k in common])
            pairs.append({
                "receptor": rec_id,
                "ligand_from": lig_id,
                "pocket_ca_rmsd": round(rmsd, 3),
            })

    # Output
    csv_path = os.path.join(OUTPUT_DIR, "pocket_pairwise_rmsd.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["receptor", "ligand_from", "pocket_ca_rmsd"])
        w.writeheader()
        w.writerows(pairs)
    print(f"Saved: {csv_path} ({len(pairs)} pairs)")
    return pairs, csv_path


def run_positive_control():
    pairs, csv_path = compute_pairwise_pocket_rmsd()

    # Load cross-docking results
    dock_path = os.path.join(OUTPUT_DIR, "cross_docking_results.csv")
    if not os.path.exists(dock_path):
        print("ERROR: cross_docking_results.csv not found")
        return
    docking = list(csv.DictReader(open(dock_path)))

    # Load B-factor data
    bf_path = os.path.join(OUTPUT_DIR, "bfactor_summary.csv")
    bf_by_pdb = {}
    if os.path.exists(bf_path):
        for r in csv.DictReader(open(bf_path)):
            bf_by_pdb.setdefault(r["pdb_id"], []).append(r)

    # Merge
    merged = []
    for pair in pairs:
        rec, lig = pair["receptor"], pair["ligand_from"]
        key = f"{rec}|{lig}"
        dock = next((d for d in docking
                     if d["receptor"] == rec and d["ligand_from"] == lig), None)
        if dock is None or dock.get("rmsd") in (None, ""):
            continue
        rmsd = float(dock["rmsd"])
        pocket_rmsd = pair["pocket_ca_rmsd"]
        bf_avg = None
        if rec in bf_by_pdb and bf_by_pdb[rec]:
            bf_avg = np.mean([float(r["avg_bfactor"]) for r in bf_by_pdb[rec]])
        merged.append({
            "receptor": rec, "ligand_from": lig,
            "docking_rmsd": rmsd,
            "pocket_ca_rmsd": pocket_rmsd,
            "bfactor_avg": bf_avg,
        })

    if len(merged) < 10:
        print(f"WARNING: only {len(merged)} matched pairs, need more")
        return

    # Stats
    pr = [m["pocket_ca_rmsd"] for m in merged]
    dr = [m["docking_rmsd"] for m in merged]
    br = [m["bfactor_avg"] for m in merged if m["bfactor_avg"] is not None]

    r_pocket, p_pocket = stats.spearmanr(pr, dr)
    print(f"\n{'='*60}")
    print("POSITIVE CONTROL: Pocket Cα RMSD vs Docking RMSD")
    print(f"{'='*60}")
    print(f"Spearman ρ = {r_pocket:.4f}, p = {p_pocket:.4f}")

    if len(br) == len(dr):
        r_bf, p_bf = stats.spearmanr(br, dr)
        print(f"\nB-factor (for comparison):")
        print(f"Spearman ρ = {r_bf:.4f}, p = {p_bf:.4f}")

    # Figure: side-by-side
    if HAS_MPL:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

        ax1.scatter(pr, dr, c="steelblue", alpha=0.6, s=25, edgecolors="gray", linewidth=0.3)
        ax1.set_xlabel("Pocket Cα RMSD (Å)", fontsize=12)
        ax1.set_ylabel("Cross-docking RMSD (Å)", fontsize=12)
        ax1.set_title(f"Pocket Similarity → Docking\nSpearman ρ={r_pocket:.3f}, p={p_pocket:.4f}", fontsize=12)

        if len(br) == len(dr):
            ax2.scatter(br, dr, c="darkorange", alpha=0.6, s=25, edgecolors="gray", linewidth=0.3)
            ax2.set_xlabel("Pocket Avg B-factor (Å²)", fontsize=12)
            ax2.set_ylabel("Cross-docking RMSD (Å)", fontsize=12)
            ax2.set_title(f"B-factor → Docking\nSpearman ρ={r_bf:.3f}, p={p_bf:.4f}", fontsize=12)

        fig.tight_layout()
        fig_path = os.path.join(FIGURES_DIR, "positive_control_comparison.png")
        fig.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\nSaved: {fig_path}")

    # Save CSV
    csv_out = os.path.join(OUTPUT_DIR, "positive_control_analysis.csv")
    with open(csv_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["receptor", "ligand_from",
                                           "docking_rmsd", "pocket_ca_rmsd", "bfactor_avg"])
        w.writeheader()
        w.writerows(merged)
    print(f"Saved: {csv_out}")


if __name__ == "__main__":
    run_positive_control()
