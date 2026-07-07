# B-factors and Docking — Group 6 Final Project

Do crystallographic **B-factors** at a protein's binding site predict how hard that
pocket is to **dock** into? This project tests that idea end to end: we survey a
protein family in the PDB, extract binding-site flexibility from B-factors, run
**cross-docking** with AutoDock Vina, and ask whether flexible pockets produce
worse (higher-RMSD) poses. We then use the same B-factor signal to guide
**ensemble docking**.

## Hypothesis

Binding-site residues with high relative B-factors are more flexible. Flexibility
makes a single rigid receptor a poorer docking target, so we expect:

1. **Cross-docking:** pairs whose receptor has higher binding-site B-factors show
   higher docking RMSD (worse accuracy).
2. **Ensemble docking:** choosing an ensemble of structures guided by B-factors
   (`BFIbs`, a B-factor Index for the binding site) beats random selection.

## Pipeline

The experiment is driven by [`experiment/main.py`](experiment/main.py), which
orchestrates five steps:

| Step | Script | What it does |
|------|--------|--------------|
| 1 | `step1_download_structures.py` | Survey the RCSB PDB for a target (e.g. CDK2) and download holo (protein + ligand) structures |
| 2 | `step2_extract_bfactors.py` | Find binding-site residues within 5 Å of the ligand and extract per-residue B-factors (relative to the pocket average) |
| 3 | `step3_cross_docking.py` | Prepare PDBQT inputs and cross-dock every ligand into every other receptor with Vina |
| 4 | `step4_analyze.py` | Test the correlation between receptor binding-site B-factors and docking RMSD |
| 5 | `step5_ensemble_docking.py` | Compute `BFIbs`, select ensembles (B-factor-guided vs. random vs. lowest-B-factor), and compare docking performance |

![Pipeline](experiment/pipeline_diagram.png)

## Usage

```bash
cd experiment

# One-time: survey the target and pick structures
python main.py survey CDK2        # prints candidate PDB IDs
# → paste the chosen IDs into TARGET_PDB_IDS in config.py
python main.py download

# Run the analysis
python main.py all                # steps 2–5 in sequence
# ...or run individual steps:
python main.py bfactors
python main.py prepare
python main.py dock
python main.py analyze
python main.py ensemble
```

Results are written to `experiment/output/` (CSVs) and `experiment/figures/` (plots).

## Prerequisites

```bash
pip install biopython numpy scipy matplotlib
```

External tools (must be on your `PATH`):
- **AutoDock Vina** — https://vina.scripps.edu/download/
- **Open Babel** — https://openbabel.org/ (PDB → PDBQT conversion)

## Configuration

All knobs live in [`experiment/config.py`](experiment/config.py): the target PDB
IDs, binding-site distance cutoff, Vina search box / exhaustiveness, the RMSD
success threshold, and the ensemble sizes and selection strategies.

## Repository layout

```
experiment/                              # the pipeline (steps 1–5 + config + diagram)
Literature_Review_B-factors_and_Docking.docx
Paper_Outline_B-factors_and_Docking.docx
```

Downloaded structures, docking results, and generated figures are produced at
runtime and are not tracked in git (see `.gitignore`).

## colab
- emily
- lola
- lucas