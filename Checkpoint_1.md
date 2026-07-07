**Checkpoint 1**

*Group 6 --- Computational Biology Final Project*

*Can B-factors identify flexible binding-site residues that reduce rigid
docking accuracy?*

# 1. Topic Introduction

Protein flexibility poses a fundamental challenge for rigid docking
algorithms: a single crystal structure captures only one snapshot of a
protein's conformational ensemble, yet different ligands may bind to
distinct conformations. Cross-docking experiments consistently reveal
that success rates drop sharply when the receptor conformation differs
from the one co-crystallized with the ligand. This project investigates
whether crystallographic B-factors---readily available in every PDB
structure---can predict which binding-site residues are flexible enough
to degrade rigid docking accuracy. We further test whether B-factors can
guide the selection of receptor ensembles for ensemble docking, using
the B-factor Index for the binding site (BFIbs) to choose which
structures to include. By comparing B-factor-guided ensemble selection
against random and lowest-B-factor baselines, we evaluate whether this
simple, zero-cost metric can improve pose prediction without the
computational expense of molecular dynamics simulations.

- **PDB Survey & B-factor Extraction** --- Download all holo structures
  of a target protein (e.g., CDK2) from the PDB; extract per-residue
  B-factors for binding-site residues within 5 Å of the co-crystallized
  ligand; compute the B-factor Index for the binding site (BFIbs) for
  each structure.

- **Cross-Docking with AutoDock Vina** --- For every pair of structures
  (receptor A + ligand B, A ≠ B), run rigid docking; measure RMSD
  between the docked pose and the crystal pose; classify each pair as
  success (RMSD \< 2 Å) or failure.

- **B-factor--Guided Ensemble Docking** --- Select receptor ensembles
  using three strategies: BFIbs-guided (structures closest to BFIbs =
  1), random (baseline), and lowest-B-factor (most rigid binding sites).
  Dock each ligand against all ensemble members, taking the best pose.
  Compare success rates to test whether B-factors improve ensemble
  selection over random.

- **Ensemble Optimizer Comparison** --- Implement an iterative ensemble
  optimizer (De Paris et al., 2018; Bolstad & Anderson, 2009) as a
  literature-based upper bound: cluster receptor conformations by
  docking-score feedback, then compare its success rate against
  B-factor-guided and random selection to benchmark our method.

# 2. Key References

1.  **Amaro, R. E., Baudry, J. L., Chodera, J. D., Demir, Ö.,
    McCammon, J. A., Miao, Y., & Smith, J. C. (2018). Ensemble docking
    in drug discovery.** *Biophysical Journal, 114*(10), 2271--2278.
    https://doi.org/10.1016/j.bpj.2018.02.038

> The definitive review of ensemble docking methodology. Describes the
> Relaxed Complex Scheme (RCS)---generating conformational ensembles
> from MD trajectories---and the core challenge: ensemble selection is
> the critical bottleneck, not ensemble generation.

2.  **Halip, L., Avram, S., & Neanu, C. (2021). The B-factor index for
    the binding site (BFIbs) to prioritize crystal protein structures
    for docking.** *Structural Chemistry, 32*(4), 1693--1699.
    https://doi.org/10.1007/s11224-021-01720-4

> Introduces BFIbs = median(B_pocket) / median(B_protein). Across 26,019
> protein--ligand complexes, BFIbs was the only crystallographic metric
> significantly correlated with docking RMSD. Structures with BFIbs \< 1
> produced the best docking results (RMSD \< 2 Å). This is the direct
> methodological predecessor to our ensemble selection strategy.

3.  **Korb, O., Olsson, T. S. G., Bowden, S. J., Hall, R. J.,
    Verdonk, M. L., Liebeschuetz, J. W., & Cole, J. C. (2012). Potential
    and limitations of ensemble docking.** *Journal of Chemical
    Information and Modeling, 52*(5), 1292--1304.
    https://doi.org/10.1021/ci300064d

> Large-scale systematic evaluation: 8 targets, 500,000 receptor
> combinations. Key finding: most randomly assembled ensembles do not
> improve enrichment over single-structure docking. Ensemble selection
> matters more than ensemble size---validating our hypothesis that a
> principled selection strategy (B-factor) should outperform random.

4.  **Palacio-Rodríguez, K., Lans, I., & Cavasotto, C. N. (2019).
    Exponential consensus ranking improves the outcome in docking and
    receptor ensemble docking.** *Scientific Reports, 9*, 5142.
    https://doi.org/10.1038/s41598-019-41594-3

> Proposes Exponential Consensus Ranking (ECR) for integrating
> multi-structure docking results. Demonstrates that ensemble docking +
> consensus ranking consistently outperforms single-structure docking
> across 4 diverse targets (CDK2, ESR1, CAH2, ADRB2). Provides a
> principled framework for how to combine results once an ensemble is
> selected.

5.  **Bottegoni, G., Kufareva, I., Totrov, M., & Abagyan, R. (2009).
    Four-dimensional docking: A fast and accurate account of discrete
    receptor flexibility in ligand docking.** *Journal of Medicinal
    Chemistry, 52*(2), 397--406. https://doi.org/10.1021/jm8009958

> Treats receptor conformations as a fourth dimension, achieving \~4×
> speedup over traditional ensemble docking while maintaining accuracy
> on 99 proteins and 300 ligands. Demonstrates that the sampling
> algorithm, not just ensemble composition, governs success.

6.  **Rao, S., Sanschagrin, P. C., Greenwood, J. R., Repasky, M. P.,
    Sherman, W., & Farid, R. (2008). Improving database enrichment
    through ensemble docking.** *Journal of Computer-Aided Molecular
    Design, 22*(9), 657--663. https://doi.org/10.1007/s10822-008-9163-6

> Demonstrated that 2--3 well-chosen receptor structures in an ensemble
> can match or exceed performance of larger ensembles. Proposed an
> ensemble selection strategy based on average GlideScore of top-ranked
> ligands, without requiring known actives.

7.  **Najmanovich, R., Kuttner, J., Sobolev, V., & Edelman, M. (2000).
    Side-chain flexibility in proteins upon ligand binding.** *Proteins:
    Structure, Function, and Bioinformatics, 39*(3), 261--268.
    https://doi.org/10.1002/(SICI)1097-0134(20000515)39:3\<261::AID-PROT90\>3.0.CO;2-4

> Found that backbone motions and side-chain rotations are decoupled---a
> fundamental limitation for using B-factors (which primarily reflect
> backbone motion) to predict side-chain flexibility. This is the most
> important counter-argument to our hypothesis.

8.  **De Paris, R., Vahl Quevedo, C., Ruiz, D. D., Gargano, F., & de
    Souza, O. N. (2018). A selective method for optimizing ensemble
    docking-based experiments on an InhA fully-flexible receptor
    model.** *BMC Bioinformatics, 19*, 235.
    https://doi.org/10.1186/s12859-018-2228-9

> Proposes a hierarchical clustering strategy for ensemble selection
> using docking-score feedback---an iterative optimizer that picks
> conformations based on how well they dock, rather than structural
> similarity alone. Serves as the literature benchmark for ensemble
> optimization.

9.  **Bolstad, A. C., & Anderson, A. C. (2009). In pursuit of virtual
    lead optimization: Pruning ensembles of receptor structures for
    increased efficiency and accuracy during docking.** *Journal of
    Computer-Aided Molecular Design, 23*(11), 755--763.
    https://doi.org/10.1007/s10822-009-9220-9

> Introduces a relative-distance-based pruning method that retains
> ensemble members whose active-site core geometry best preserves the
> conserved ligand-binding pharmacophore. Demonstrates that geometric
> pruning can match the accuracy of full ensembles at a fraction of the
> computational cost.

# 3. Major Blocks & Challenges

The following challenges represent the most significant obstacles we
anticipate, prioritized by risk to the project's core claims.

## 3.1 Backbone--Sidechain Decoupling (Highest Risk)

Backbone B-factors do not directly report on side-chain rotations
(Najmanovich et al., 2000)---a residue with low backbone B-factor may
still have a flexible side-chain that blocks the pocket. We will
separately analyze residues where side-chain flexibility matters (e.g.,
gatekeeper residues in kinases) vs. backbone-only motion.

## 3.2 Target Selection Determines Signal-to-Noise

Too little conformational diversity produces a ceiling effect (all
docking succeeds); too much produces a floor effect (all docking fails).
CDK2 structures solved with similar inhibitors may lack meaningful
variation. We will compute pairwise RMSD across structures to validate
diversity before committing to a target.

## 3.3 Statistical Power with Small Ensembles

With limited structurally distinct conformations, BFIbs selection
operates on few meaningful options---reducing statistical power. We will
report effect sizes and confidence intervals alongside p-values,
emphasizing qualitative trends over significance alone.

## 3.4 AutoDock Vina's Limitations

Vina is faster but less accurate than Glide, GOLD, or ICM (Korb et al.,
2012); its scoring function may wash out subtle conformational effects.
Cross-validation with a second scoring function is planned if time
permits, and we acknowledge that conclusions may not generalize to all
docking engines.

## 3.5 Correlational vs. Causal Claims

We can only claim correlation, not causation. Confounds include
resolution bias, crystal packing artifacts, and ligand-induced B-factor
changes. We normalize B-factors via BFIbs, check for resolution
dependence, and explicitly distinguish correlation from causation
throughout.

# 4. Work split

  ------------------------------------------------------------------------
  Name            Percent      Contributions
  --------------- ------------ -------------------------------------------
  Lucas           40           Research, coming up with ideas

  James           40           Research, coming up with ideas

  Lola            10           Read papers

  Emily           10           Read papers
  ------------------------------------------------------------------------
