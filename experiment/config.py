"""
Configuration for the B-factor / docking experiment.
=====================================================
Change TARGET_PDB_IDS to select a different protein target.
All paths are relative to this file's directory.
"""

import os

# === Directory structure ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
VINA_DIR = os.path.join(DATA_DIR, "vina_inputs")   # PDBQT files
RESULTS_DIR = os.path.join(DATA_DIR, "docking_results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

for d in [PDB_DIR, VINA_DIR, RESULTS_DIR, FIGURES_DIR, OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# === Target selection ===
# CDK2: 49+ holo structures, well-studied, moderate flexibility
# Change this list after running step1 to survey available structures.
TARGET_PDB_IDS = [
    '5MHQ', '3RM6', '3PXZ', '3QQJ', '3R73', '3RAI',
    '3PY0', '3QQF', '3QRT', '3QWJ', '3QX2', '3QXO',
    '3QZG', '3QZI', '3R28', '3R6X', '3R71', '3R7U',
    '3R83', '3ROY', '3RPO', '3PXY', '3PXF', '3R1S',
    '3R1Y', '3R8M', '3R8P', '3QQL', '3QWK', '3R1Q',
]

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
