"""
Configuration for the B-factor / docking experiment.
=====================================================
Change TARGET_PDB_IDS to select a different protein target.
All paths are relative to this file's directory.
"""

import os
import re

# === Target selection ===
# Set BFIBS_TARGET_NAME in the shell to run another target while keeping
# generated files separated by protein name.
TARGET_NAME = os.environ.get("BFIBS_TARGET_NAME", "CDK2")
TARGET_SLUG = re.sub(r"[^A-Za-z0-9_.-]+", "_", TARGET_NAME).strip("_")

# === Directory structure ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data", TARGET_SLUG)
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
VINA_DIR = os.path.join(DATA_DIR, "vina_inputs")   # PDBQT files
RESULTS_DIR = os.path.join(DATA_DIR, "docking_results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures", TARGET_SLUG)
OUTPUT_DIR = os.path.join(BASE_DIR, "output", TARGET_SLUG)

for d in [PDB_DIR, VINA_DIR, RESULTS_DIR, FIGURES_DIR, OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# === Target PDB IDs ===
TARGET_PDB_IDS_BY_NAME = {
    # CDK1: Cyclin-Dependent Kinase 1, 8 holo X-ray structures (limited availability)
    "CDK1": [
        '6GU2', '6GU3', '6GU4', '6GU6', '6GU7',
        '4Y72', '5HQ0', '5LQF',
    ],
    # CDK2: 49+ holo structures, well-studied, moderate flexibility
    "CDK2": [
        '5MHQ', '3RM6', '3PXZ', '3QQJ', '3R73', '3RAI',
        '3PY0', '3QQF', '3QRT', '3QWJ', '3QX2', '3QXO',
        '3QZG', '3QZI', '3R28', '3R6X', '3R71', '3R7U',
        '3R83', '3ROY', '3RPO', '3PXY', '3PXF', '3R1S',
        '3R1Y', '3R8M', '3R8P', '3QQL', '3QWK', '3R1Q',
    ],
    # Trypsin: top X-ray structures from the RCSB full-text survey.
    "Trypsin": [
        '1XVO', '1HJ8', '1XVM', '1PPZ', '1PQA', '2A31',
        '5PTP', '1UTM', '1TRY', '1MCT', '1A0J', '1TPO',
        '1SGT', '3PTB', '1AVW', '1F2S', '1EPT', '1AKS',
        '1MBQ', '1UTJ', '1BIT', '7JR2', '1AVX', '6DWF',
        '1YF4', '6DWH', '2FPZ', '7JR1', '1UHB', '1TX6',
    ],
}

TARGET_PDB_IDS = TARGET_PDB_IDS_BY_NAME.get(TARGET_NAME, TARGET_PDB_IDS_BY_NAME["CDK2"])

# === Binding site parameters ===
BINDING_SITE_DISTANCE = 5.0   # Å from co-crystallized ligand heavy atoms

# === Vina parameters ===
VINA_BINARY = os.path.join(BASE_DIR, "vina")  # Bundled Vina binary (auto-detected)
VINA_BOX_SIZE = (25, 25, 25)  # Å — search box dimensions
VINA_EXHAUSTIVENESS = 64      # Higher = more thorough but slower
VINA_NUM_MODES = 9            # Number of binding poses to generate

# === Analysis thresholds ===
RMSD_SUCCESS_THRESHOLD = 2.0  # Å — pose considered correct if RMSD < this
BFACTOR_HIGH_PERCENTILE = 75  # Top quartile = "high B-factor" residues

# === Ensemble docking parameters ===
ENSEMBLE_SIZES = [2, 3, 5]          # Test different ensemble sizes
ENSEMBLE_STRATEGIES = [
    "bfactor_guided",    # Select structures with BFIbs closest to 1
    "random",            # Random selection (baseline)
    "lowest_bfactor",    # Select structures with lowest BFIbs (most rigid pocket)
    "enopt",             # ML-based ensemble selection via EnOpt (Bhatt et al., 2024)
]
NUM_RANDOM_ENSEMBLES = 20   # Number of random ensembles for statistical comparison

# === Visualization ===
FIGURE_DPI = 150
FIGURE_STYLE = "seaborn-v0_8-whitegrid"
