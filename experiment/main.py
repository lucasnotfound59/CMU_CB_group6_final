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
    python main.py ensemble             # Run ensemble docking (all strategies)
    python main.py enopt                # Run EnOpt ML-based ensemble selection
    python main.py all                  # Run everything (except survey)

Prerequisites:
    pip install biopython numpy scipy matplotlib pandas
    # Install AutoDock Vina: https://vina.scripps.edu/download/
    # Install Open Babel: https://openbabel.org/ (for PDB → PDBQT conversion)
    # EnOpt deps: pip install scikit-learn xgboost plotly
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))


def cmd_survey(target):
    import re
    from step1_download_structures import survey_target

    results = survey_target(target)
    if not results:
        print("No valid structures found.")
        return

    # Build PDB ID list (top 30, sorted by resolution)
    pdb_ids = [r["pdb_id"] for r in results[:30]]

    # Auto-update config.py
    config_path = os.path.join(os.path.dirname(__file__), "config.py")
    with open(config_path, "r") as f:
        content = f.read()

    # Format the new TARGET_PDB_IDS block
    lines = ["TARGET_PDB_IDS = [\n"]
    for i in range(0, len(pdb_ids), 6):
        chunk = ", ".join(repr(pid) for pid in pdb_ids[i:i+6])
        lines.append(f"    {chunk},\n")
    lines.append("]")
    new_block = "".join(lines)

    # Replace the TARGET_PDB_IDS block in config.py
    pattern = r"TARGET_PDB_IDS = \[[^\]]*\]"
    content = re.sub(pattern, new_block, content, flags=re.DOTALL)

    with open(config_path, "w") as f:
        f.write(content)

    print(f"\n✓ Updated config.py with {len(pdb_ids)} PDB IDs")


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
    """Run ensemble docking with all selection strategies."""
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
    
    # Select ensembles (includes enopt if enopt_ensemble.csv exists)
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


def cmd_enopt():
    """Run EnOpt ML-based ensemble selection benchmark."""
    import step6_enopt
    
    print("=" * 60)
    print("ENOPT BENCHMARK")
    print("=" * 60)
    
    # Export matrix
    print("\nStep 1: Export EnOpt score matrix...")
    matrix_path = step6_enopt.export_enopt_matrix()
    known_path = step6_enopt.create_known_ligands(matrix_path)
    
    # Run EnOpt
    print("\nStep 2: Running EnOpt...")
    out_prefix = step6_enopt.run_enopt(matrix_path, known_path)
    
    # Parse ensemble
    print("\nStep 3: Parsing EnOpt ensemble selection...")
    step6_enopt.parse_enopt_ensemble(out_prefix)
    
    print("\n" + "=" * 60)
    print("ENOPT BENCHMARK COMPLETE")
    print("=" * 60)
    print("\nRun 'python main.py ensemble' to include EnOpt in strategy comparison.")


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
    print("STEP 6: EnOpt ML benchmark")
    print("=" * 60)
    cmd_enopt()

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
    elif command == "enopt":
        cmd_enopt()
    elif command == "all":
        cmd_all()
    else:
        print(f"Unknown command: {command}")
        print(__doc__)
