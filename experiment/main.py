#!/usr/bin/env python3
"""
Main pipeline: B-factor → Docking Accuracy Experiment
======================================================
Orchestrates all steps of the computational experiment.

Usage:
    python main.py survey CDK2          # Find CDK2 structures in PDB
    python main.py download             # Download structures listed in config.py
    python main.py bfactors             # Extract binding-site B-factors
    python main.py prepare              # Prepare Vina inputs (PDBQT files)
    python main.py dock                 # Run cross-docking
    python main.py analyze              # Analyze B-factor vs RMSD correlation
    python main.py ensemble             # Run ensemble docking (B-factor guided)
    python main.py all                  # Run everything (except survey)

Prerequisites:
    pip install biopython numpy scipy matplotlib
    # Install AutoDock Vina: https://vina.scripps.edu/download/
    # Install Open Babel: https://openbabel.org/ (for PDB → PDBQT conversion)
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))


def cmd_survey(target):
    from step1_download_structures import survey_target
    survey_target(target)


def cmd_download():
    from step1_download_structures import download_structures
    from config import TARGET_PDB_IDS
    if not TARGET_PDB_IDS:
        print("ERROR: No PDB IDs in config.TARGET_PDB_IDS")
        print("Run: python main.py survey <target_name>  first")
        sys.exit(1)
    download_structures(TARGET_PDB_IDS)


def cmd_bfactors():
    import step2_extract_bfactors
    step2_extract_bfactors.main()


def cmd_prepare():
    from step3_cross_docking import prepare_all
    prepare_all()


def cmd_dock():
    from step3_cross_docking import run_cross_docking
    run_cross_docking()


def cmd_analyze():
    import step4_analyze
    step4_analyze.analyze()


def cmd_ensemble():
    """Run ensemble docking with B-factor guided selection."""
    import step5_ensemble_docking as ens
    
    print("=" * 60)
    print("ENSEMBLE DOCKING PIPELINE")
    print("=" * 60)
    
    # Calculate BFIbs
    print("\nStep 1: Calculating BFIbs scores...")
    bfibs_data = ens.compute_all_bfibs()
    if not bfibs_data:
        print("ERROR: No BFIbs data. Run bfactors step first.")
        return
    
    # Select ensembles
    print("\nStep 2: Selecting ensembles by different strategies...")
    ensembles = ens.select_all_ensembles(bfibs_data)
    ens.save_ensemble_selections(ensembles, bfibs_data)
    
    # Run ensemble docking
    print("\nStep 3: Running ensemble docking...")
    results = ens.run_ensemble_docking(ensembles)
    
    # Compare strategies
    print("\nStep 4: Comparing strategies...")
    ens.compare_strategies(results)
    
    print("\n" + "=" * 60)
    print("ENSEMBLE PIPELINE COMPLETE")
    print("=" * 60)


def cmd_all():
    print("=" * 60)
    print("STEP 1: Extract B-factors")
    print("=" * 60)
    cmd_bfactors()

    print("\n" + "=" * 60)
    print("STEP 2: Prepare Vina inputs")
    print("=" * 60)
    cmd_prepare()

    print("\n" + "=" * 60)
    print("STEP 3: Run cross-docking")
    print("=" * 60)
    cmd_dock()

    print("\n" + "=" * 60)
    print("STEP 4: Analyze results")
    print("=" * 60)
    cmd_analyze()

    print("\n" + "=" * 60)
    print("STEP 5: Ensemble docking")
    print("=" * 60)
    cmd_ensemble()

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print("\nCheck output/ for CSV results")
    print("Check figures/ for analysis plots")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    command = sys.argv[1].lower()

    if command == "survey":
        if len(sys.argv) < 3:
            print("Usage: python main.py survey <target_name>")
            sys.exit(1)
        cmd_survey(sys.argv[2])
    elif command == "download":
        cmd_download()
    elif command == "bfactors":
        cmd_bfactors()
    elif command == "prepare":
        cmd_prepare()
    elif command == "dock":
        cmd_dock()
    elif command == "analyze":
        cmd_analyze()
    elif command == "ensemble":
        cmd_ensemble()
    elif command == "all":
        cmd_all()
    else:
        print(f"Unknown command: {command}")
        print(__doc__)
