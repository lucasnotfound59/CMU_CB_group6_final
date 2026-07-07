**Project Plan**

*Group 6 --- Computational Biology Final Project*

Research Question

Can B-factors identify flexible binding-site residues that reduce rigid
docking accuracy? We also test whether B-factors can guide ensemble
docking selection via the B-factor Index for the binding site (BFIbs).

Pipeline

Five steps, driven by experiment/main.py

- **PDB Survey & B-factor Extraction ---** Download holo structures of a
  target (e.g., CDK2); extract per-residue B-factors within 5Å of
  ligand; compute BFIbs.

- **Cross-Docking with Vina ---** For every receptor A + ligand B (A≠B),
  run rigid docking; measure RMSD vs crystal pose; classify success
  (RMSD\<2Å) or failure.

- **Correlation Analysis ---** Spearman ρ between receptor B-factors and
  docking RMSD; t-test comparing B-factors of successful vs failed
  pairs.

- **B-factor--Guided Ensemble Docking ---** Select ensembles by BFIbs
  (closest to 1), random (baseline), and lowest BFIbs (most rigid).
  Compare success rates.

- **Ensemble Optimizer Benchmark ---** Implement literature optimizer
  (De Paris 2018; Bolstad 2009) as an upper bound for comparison.

Deliverables

- Paper (6 sections: Intro, Background, Evidence, Experiment,
  Discussion, Conclusion)

- Presentation slides

- Experiment code + results (CSV in output/, figures in figures/)

Timeline (14 Days)

- **Day 1--2:** Install Vina/obabel; survey PDB; select target; download
  structures

- **Day 3--4:** Extract B-factors; compute BFIbs; prepare PDBQT files

- **Day 5--7:** Run cross-docking + ensemble docking; statistical
  analysis; generate figures

- **Day 8--10:** Write paper (Intro through Conclusion)

- **Day 11--12:** Create presentation slides

- **Day 13--14:** Revise paper; rehearse presentation

Key References

- Amaro et al. (2018) --- Ensemble docking in drug discovery. Biophys.
  J.

- Halip et al. (2021) --- BFIbs to prioritize structures for docking.
  Struct. Chem.

- Korb et al. (2012) --- Potential and limitations of ensemble docking.
  JCIM.

- Palacio-Rodríguez et al. (2019) --- ECR improves ensemble docking.
  Sci. Rep.

- Bottegoni et al. (2009) --- 4D docking. J. Med. Chem.

- Rao et al. (2008) --- Improving enrichment through ensemble docking.
  JCAMD.

- Najmanovich et al. (2000) --- Side-chain flexibility. Proteins.

- De Paris et al. (2018) --- Selective ensemble optimization. BMC
  Bioinf.

- Bolstad & Anderson (2009) --- Pruning ensembles. JCAMD.

Key Risks

- **Backbone--sidechain decoupling:** B-factors measure backbone motion,
  not side-chain rotations (Najmanovich 2000).

- **Target selection:** Need structural diversity---too little gives
  ceiling effects, too much gives floor effects.

- **Statistical power:** Limited distinct conformations may reduce
  power. Report effect sizes + CIs.

- **Vina limitations:** Less accurate than Glide/GOLD/ICM.
  Cross-validate if time permits.

- **Correlation ≠ causation:** Confounds: resolution bias, crystal
  packing, ligand-induced B-factor changes.
