#!/usr/bin/env python3
"""
Step 1: Survey and download protein structures from the PDB.
=============================================================
Usage:
    python step1_download_structures.py --survey CDK2     # Find all CDK2 structures
    python step1_download_structures.py --download         # Download structures in config

This script:
    1. Queries the RCSB PDB API for structures of a given target
    2. Filters for holo structures (protein + ligand)
    3. Downloads PDB files for cross-docking analysis
"""

import argparse
import json
import os
import ssl
import sys
import urllib.request
import urllib.parse

import certifi
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

sys.path.insert(0, os.path.dirname(__file__))
from config import PDB_DIR, TARGET_PDB_IDS

HTTPS_CONTEXT = ssl.create_default_context(cafile=certifi.where())


# === RCSB API ===

def search_pdb(query: str, rows: int = 500) -> list[dict]:
    """
    Query RCSB PDB search API.
    Returns list of {pdb_id, title, resolution, experimental_method}.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query}
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}}
    }
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, context=HTTPS_CONTEXT) as resp:
        data = json.loads(resp.read())

    results = []
    for hit in data.get("result_set", []):
        pdb_id = hit["identifier"].upper()
        results.append({"pdb_id": pdb_id})
    return results


def get_structure_info(pdb_id: str) -> dict:
    """
    Fetch structure metadata from RCSB GraphQL API.
    Returns resolution, title, experimental method.
    """
    url = "https://data.rcsb.org/graphql"
    query = f"""
    {{
        entry(entry_id: "{pdb_id}") {{
            rcsb_entry_info {{
                resolution_combined
                experimental_method
            }}
            struct {{
                title
            }}
        }}
    }}
    """
    payload = {"query": query}
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, context=HTTPS_CONTEXT) as resp:
        data = json.loads(resp.read())

    entry = data["data"]["entry"]
    info = entry["rcsb_entry_info"]
    return {
        "pdb_id": pdb_id,
        "title": entry["struct"]["title"],
        "resolution": info["resolution_combined"],
        "method": info["experimental_method"],
    }


# === PDB Download ===

def download_pdb(pdb_id: str, output_dir: str) -> str:
    """
    Download a PDB file from RCSB. Returns the file path.
    """
    pdb_id_lower = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id_lower}.pdb"
    out_path = os.path.join(output_dir, f"{pdb_id_lower}.pdb")

    if os.path.exists(out_path):
        print(f"  [skip] {pdb_id} already downloaded")
        return out_path

    print(f"  [download] {pdb_id} ...")
    with urllib.request.urlopen(url, context=HTTPS_CONTEXT) as resp:
        with open(out_path, "wb") as f:
            f.write(resp.read())
    return out_path


# === Binding Site Extraction ===

class LigandSelect(Select):
    """Select only protein chains + a specific ligand."""
    def __init__(self, ligand_resname: str):
        self.ligand_resname = ligand_resname
        self.protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }

    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        if resname in self.protein_residues:
            return True
        if resname == self.ligand_resname:
            return True
        return False


def extract_binding_site_residues(pdb_path: str, ligand_resname: str,
                                   distance_cutoff: float = 5.0) -> list:
    """
    Parse a PDB file and return binding-site residues within `distance_cutoff`
    Å of the ligand.

    Returns:
        list of dicts: [{chain, resseq, resname, bfactor, ca_bfactor}, ...]
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    model = structure[0]

    # Find ligand atoms
    ligand_atoms = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip() == ligand_resname:
                for atom in residue:
                    ligand_atoms.append(atom.get_coord())

    if not ligand_atoms:
        return []

    import numpy as np
    ligand_coords = np.array(ligand_atoms)

    # Find protein residues near the ligand
    binding_site = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname not in LigandSelect(ligand_resname).protein_residues:
                continue

            # Check if any atom is within distance cutoff
            is_near = False
            for atom in residue:
                dists = np.linalg.norm(ligand_coords - atom.get_coord(), axis=1)
                if np.min(dists) < distance_cutoff:
                    is_near = True
                    break

            if is_near:
                # Extract B-factors
                atom_bfactors = [atom.get_bfactor() for atom in residue
                                 if atom.get_id() == "CA"]
                ca_bfactor = atom_bfactors[0] if atom_bfactors else None
                avg_bfactor = sum(a.get_bfactor() for a in residue) / len(residue)

                binding_site.append({
                    "chain": chain.get_id(),
                    "resseq": residue.get_id()[1],
                    "resname": resname,
                    "avg_bfactor": round(avg_bfactor, 2),
                    "ca_bfactor": round(ca_bfactor, 2) if ca_bfactor else None,
                })

    return binding_site


# === Main ===

def survey_target(target_name: str):
    """Survey available structures for a target protein."""
    print(f"\n{'='*60}")
    print(f"Surveying PDB for: {target_name}")
    print(f"{'='*60}\n")

    # Search
    results = search_pdb(target_name)
    print(f"Found {len(results)} total entries\n")

    # Get info for each (limit to first 50 for speed)
    infos = []
    for item in results[:50]:
        try:
            info = get_structure_info(item["pdb_id"])
            if info["resolution"] and info["resolution"][0] is not None:
                infos.append(info)
        except Exception as e:
            print(f"  [warn] Could not fetch info for {item['pdb_id']}: {e}")

    # Filter: X-ray, resolution <= 3.0Å
    def valid(info):
        if not info.get("method"):
            return False
        if "X-RAY" not in str(info["method"]).upper():
            return False
        res = info.get("resolution")
        if res is None:
            return False
        res_val = res[0] if isinstance(res, list) else res
        return res_val is not None and res_val <= 3.0

    filtered = [i for i in infos if valid(i)]

    print(f"After filtering (X-ray, ≤3.0Å): {len(filtered)} structures\n")

    # Sort by resolution
    filtered.sort(key=lambda x: x["resolution"][0])

    # Print table
    print(f"{'PDB ID':<10} {'Resolution':<12} {'Title'}")
    print("-" * 70)
    for info in filtered[:30]:
        res_val = info['resolution'][0] if isinstance(info['resolution'], list) else info['resolution']
        res = f"{res_val:.2f} Å" if res_val else "N/A"
        title = info['title'][:45]
        print(f"{info['pdb_id']:<10} {res:<12} {title}")

    # Save list for config.py
    pdb_ids = [i["pdb_id"] for i in filtered[:30]]
    print(f"\n\n# Add these to config.py TARGET_PDB_IDS:")
    print(f"TARGET_PDB_IDS = {pdb_ids}")

    return filtered


def download_structures(pdb_ids: list[str]):
    """Download PDB files for the given IDs."""
    print(f"\n{'='*60}")
    print(f"Downloading {len(pdb_ids)} structures")
    print(f"{'='*60}\n")

    downloaded = []
    for pdb_id in pdb_ids:
        try:
            path = download_pdb(pdb_id, PDB_DIR)
            downloaded.append({"pdb_id": pdb_id, "path": path})
        except Exception as e:
            print(f"  [error] {pdb_id}: {e}")

    print(f"\nDownloaded {len(downloaded)}/{len(pdb_ids)} structures")
    return downloaded


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download PDB structures")
    parser.add_argument("--survey", type=str, help="Survey structures for a target name")
    parser.add_argument("--download", action="store_true", help="Download structures in config")
    args = parser.parse_args()

    if args.survey:
        survey_target(args.survey)
    elif args.download:
        if not TARGET_PDB_IDS:
            print("ERROR: No PDB IDs in config.TARGET_PDB_IDS. Run --survey first.")
            sys.exit(1)
        download_structures(TARGET_PDB_IDS)
    else:
        parser.print_help()
