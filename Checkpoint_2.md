**Checkpoint 2 — Code Check-in**

*Group 6 — Computational Biology Final Project*

*Can B-factors identify flexible binding-site residues that reduce rigid docking accuracy?*

**Draft Code:** [github.com/lucasnotfound59/CMU_CB_group6_final](https://github.com/lucasnotfound59/CMU_CB_group6_final) (`experiment/` directory)

# 1. Model Application: Three-Step Pipeline

We extend the BFIbs framework (Halip et al., 2021) from binding-site-level
to residue-level prediction using publicly available PDB structures and
AutoDock Vina. The pipeline is organized as five coordinated scripts
orchestrated by `experiment/main.py`, grouped into three logical phases:

## Phase 1: Data Acquisition & B-factor Extraction

**Scripts:** `step1_download_structures.py`, `step2_extract_bfactors.py`

- Query the [RCSB PDB Search API](https://search.rcsb.org) for all holo
  structures of a target protein (e.g., CDK2) using `urllib` — no API key
  required, data is fully public.
- Filter by experimental method (X-ray, ≤3.0 Å resolution) and download
  `.pdb` files via `urllib.urlretrieve`.
- Identify the co-crystallized ligand heuristically (largest non-standard,
  non-solvent residue with ≥5 heavy atoms; excludes HOH, SO₄, GOL, etc.).
- Compute per-residue B-factors for all protein residues within 5 Å of the
  ligand using BioPython's `PDBParser`.
- Calculate **relative B-factor** = B(residue) / B(pocket average), analogous
  to BFIbs (Halip et al., 2021) but at the residue level.
- Output: `output/bfactor_summary.csv` (per-residue), `output/bfactor_per_structure.csv`.

**Helper functions:**
| Function | Role |
|----------|------|
| `search_pdb(query)` | RCSB Search API v2 — returns all matching PDB IDs |
| `get_structure_info(pdb_id)` | RCSB GraphQL API — fetches resolution, method, title |
| `download_pdb(pdb_id)` | Downloads `.pdb` file, skips if cached |
| `identify_ligand(pdb_path)` | Heuristic ligand detection (largest non-AA3, non-solvent) |
| `extract_bfactors(pdb_path, ligand, cutoff)` | Returns list of binding-site residue B-factor dicts |
| `compute_relative_bfactors(residues)` | Normalizes per-residue B-factors by pocket mean |

## Phase 2: Cross-Docking with AutoDock Vina

**Script:** `step3_cross_docking.py`

- For every pair of structures (receptor *i*, ligand *j*, *i* ≠ *j*):
  1. **Prepare receptor**: strip all non-protein atoms from PDB *i* via
     BioPython `PDBIO` with a `ProteinOnly` selector.
  2. **Prepare ligand**: extract only ligand *j* atoms, writing a minimal
     HETATM PDB with proper ATOM serial numbering.
  3. **Convert PDB → PDBQT**: use Open Babel (`obabel`) or MGLTools
     (`prepare_receptor4.py` / `prepare_ligand4.py`). Gasteiger charges are
     assigned for ligands; receptor is protonated with `-xr -h`.
  4. **Define search box**: geometric center of the co-crystallized ligand +
     configurable box size (default 25×25×25 Å³).
  5. **Run Vina**: `subprocess.run` with exhaustiveness=64, num_modes=9.
  6. **Compute RMSD**: heavy-atom RMSD between best docked pose and the
     crystal pose of ligand *j* (BioPython `PDBParser`, `numpy`).
- Output: `output/cross_docking_results.csv` (receptor, ligand_from, rmsd, affinity).

**Helper functions:**
| Function | Role |
|----------|------|
| `prepare_receptor_pdb(pdb, ligand, out)` | Strips non-protein atoms via BioPython selector |
| `prepare_ligand_pdb(pdb, ligand, out)` | Extracts ligand HETATM records with correct formatting |
| `pdb_to_pdbqt(pdb, pdbqt, is_receptor)` | Converts via `obabel` (fallback: MGLTools scripts) |
| `get_binding_site_center(pdb, ligand)` | Returns (x, y, z) of ligand geometric center |
| `identify_ligand(pdb)` | Same heuristic as step2 for consistency |
| `run_vina(rec_pdbqt, lig_pdbqt, out, center, box)` | Subprocess call to Vina, parses affinity |
| `compute_rmsd(ref_pdb, docked_pdb, ligand)` | Heavy-atom RMSD between reference and docked pose |

## Phase 3: Ensemble Docking & Strategy Comparison

**Scripts:** `step4_analyze.py`, `step5_ensemble_docking.py`

- **Correlation analysis (step 4):** Merge B-factor data with cross-docking
  results. For each docking pair, compute receptor pocket statistics (mean,
  max, std B-factor). Run Pearson/Spearman correlation and a t-test comparing
  B-factors of successful (RMSD < 2 Å) vs. failed pairs. Generate
  scatter plots and distribution histograms via matplotlib.
- **BFIbs computation (step 5):** BFIbs = median(B_pocket) / median(B_protein)
  for each structure (Halip et al., 2021).
- **Ensemble selection:** Three strategies — (a) BFIbs-guided: structures
  with BFIbs closest to 1.0; (b) random: 20 ensembles for statistical
  comparison; (c) lowest-B-factor: most rigid pockets.
- **Ensemble evaluation:** For each ligand, dock against all ensemble
  members and take the best RMSD. Compare success rates across strategies
  at ensemble sizes n = {2, 3, 5}.
- **Benchmark comparison:** The Ensemble Optimizer framework (Bhatt et al.,
  2024; De Paris et al., 2018; Bolstad & Anderson, 2009) defines the
  literature upper bound — these use docking-score feedback or geometric
  pruning instead of crystallographic B-factors.

**Helper functions (step 5):**
| Function | Role |
|----------|------|
| `calculate_bfibs(pdb_id)` | Computes BFIbs = median(pocket_B) / median(protein_B) |
| `compute_all_bfibs()` | BFIbs for all structures → `output/bfibs_scores.csv` |
| `select_ensemble_bfactor_guided(data, size)` | Selects structures with BFIbs closest to 1.0 |
| `select_ensemble_lowest_bfactor(data, size)` | Selects most rigid binding sites |
| `select_ensemble_random(data, size)` | Random baseline selection |
| `select_all_ensembles(data)` | Generates ensembles for all strategies × sizes |
| `run_ensemble_docking(ensembles)` | Evaluates ensemble performance using cross-docking results |
| `compare_strategies(results)` | Statistical comparison + bar chart via matplotlib |

## Configuration & Organization

All parameters live in `experiment/config.py`: target PDB IDs, binding-site
distance cutoff (5 Å), Vina box size / exhaustiveness / num_modes, RMSD
success threshold (2 Å), ensemble sizes [2, 3, 5], and number of random
ensembles (20). Directory structure is auto-created:

```
experiment/
├── main.py                       # CLI orchestrator
├── config.py                     # All tunable parameters
├── step1_download_structures.py  # PDB survey + download
├── step2_extract_bfactors.py     # Per-residue B-factor extraction
├── step3_cross_docking.py        # Vina cross-docking + RMSD
├── step4_analyze.py              # Statistical correlation analysis
├── step5_ensemble_docking.py     # BFIbs + ensemble selection + comparison
├── data/pdb_files/               # Downloaded .pdb structures
├── data/vina_inputs/             # PDBQT receptor/ligand files
├── data/docking_results/         # Vina output PDBQT
├── output/                       # All CSVs
└── figures/                      # All plots
```

# 2. Major Blocks & Challenges (Code/Application)

## 2.1 Vina Binary Dependency (Resolved)

AutoDock Vina cannot be installed via `pip` — it requires a pre-compiled
binary. The pipeline checks for `vina` on `PATH` at runtime; `config.py`
allows specifying a custom binary path (`VINA_BINARY`). Users without Vina
installed cannot run steps 3–5.

## 2.2 PDB → PDBQT Conversion Reliability (Active)

Open Babel (`obabel`) is the primary converter, with MGLTools scripts as
fallback. However, PDBQT conversion can fail silently for non-standard
ligand chemistries (e.g., boron-containing inhibitors, covalent adducts).
We log failures but do not yet validate PDBQT output automatically.

## 2.3 Ligand Identification Heuristics

Both step2 and step3 independently identify the co-crystallized ligand via
the same heuristic (largest non-AA3, non-solvent residue). This works for
most kinase inhibitors but fails for: (a) multi-ligand complexes (cofactors
like ATP + drug), (b) peptide ligands misidentified as protein, (c) PDB
entries where the "drug" is smaller than a crystallization additive.
We currently skip ambiguous cases rather than guessing.

## 2.4 Statistical Power with Small Ensembles

With N structures (likely 10–15 after filtering), the BFIbs-guided ensemble
selects only one combination per size. Random ensembles use 20 draws for
statistical comparison, but the guided strategy has n=1. This limits
statistical inference to qualitative trends — we report effect sizes and
confidence intervals rather than relying solely on p-values.

## 2.5 RMSD Computation Correctness

`compute_rmsd()` aligns the first *n* atoms without Kabsch superposition,
relying on the assumption that Vina output atoms are in the same order as
the reference PDB. This is fragile — if Vina reorders atoms or the
reference and docked ligand have different heavy-atom counts, the RMSD is
computed on a truncated atom set. We plan to add RDKit-based RMSD with
maximum common substructure alignment if time permits.

## 2.6 Cross-Docking Is O(N²)

For N structures, the number of cross-docking pairs is N × (N−1). At ~60
seconds per Vina run (exhaustiveness=64), 12 structures produce 132 pairs
≈ 2.2 hours of wall-clock time. This is manageable for CDK2 but may become
a bottleneck for larger targets or higher exhaustiveness.

# 3. Work Split

| Name | Percent | Contributions |
|------|---------|---------------|
| Lucas | 40 | Pipeline architecture, step3 cross-docking, step5 ensemble docking |
| James | 40 | Step1 PDB survey, step2 B-factor extraction, step4 analysis |
| Lola | 10 | Literature review, paper outline, reference formatting |
| Emily | 10 | Test data preparation, presentation slides |

# References

1. Halip, L., Avram, S., & Neanu, C. (2021). The B-factor index for the binding site (BFIbs). *Structural Chemistry*, 32(4), 1693–1699.
2. Bhatt, R., Wang, A., & Durrant, J. D. (2024). Teaching old docks new tricks with machine learning enhanced ensemble docking. *Scientific Reports*, 14, 22489.
3. De Paris, R., et al. (2018). A selective method for optimizing ensemble docking-based experiments. *BMC Bioinformatics*, 19, 235.
4. Bolstad, A. C., & Anderson, A. C. (2009). In pursuit of virtual lead optimization: Pruning ensembles. *Journal of Computer-Aided Molecular Design*, 23(11), 755–763.
5. Amaro, R. E., et al. (2018). Ensemble docking in drug discovery. *Biophysical Journal*, 114(10), 2271–2278.
6. Korb, O., et al. (2012). Potential and limitations of ensemble docking. *Journal of Chemical Information and Modeling*, 52(5), 1292–1304.
