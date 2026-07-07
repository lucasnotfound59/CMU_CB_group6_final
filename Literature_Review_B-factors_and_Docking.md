*Can B-factors identify flexible binding-site residues\
that reduce rigid docking accuracy?*

# Research Question Overview

This document compiles the key literature supporting the research
question: Can B-factors (crystallographic temperature factors) identify
flexible binding-site residues that reduce the accuracy of rigid-body
docking predictions? The literature is organized into four thematic
sections, each with annotated entries including citation, summary,
applicability to our research, and open questions.

# Section 1: Protein Flexibility in Docking --- Comprehensive Reviews

## 1.1 Lexa & Carlson (2012)

**Citation:** Lexa, K. W., & Carlson, H. A. (2012). Protein flexibility
in docking and surface mapping. Quarterly Reviews of Biophysics, 45(3),
301--343. https://doi.org/10.1017/S0033583512000066

**Summary / Main Point:** A comprehensive review (179 references)
systematically cataloging how protein flexibility affects docking and
surface mapping. The paper classifies methods for incorporating
flexibility into four categories: (1) Refinement (post-docking
adjustment), (2) Average/unified structures, (3) Serial docking to
ensemble structures, and (4) On-the-fly conformational changes during
docking. Extensive tables document specific methods (soft docking,
rotamer libraries, induced fit docking, ensemble docking) with
quantitative RMSD results across multiple target proteins. Key finding:
Nichols et al. (2011) demonstrated that the predictive power of
MD-derived conformations is negatively correlated with binding site
flexibility --- highly flexible sites are less useful for docking.

**What We Can Use:** • Provides the complete methodological landscape
for handling receptor flexibility in docking\
• Tables 1--6 compile success/failure data with RMSD values for each
method across targets\
• The cross-docking problem is clearly illustrated (Figures 1--2): same
protein + different ligands → different pocket shapes\
• Directly supports the argument that identifying flexible residues a
priori is the unsolved bottleneck

**What to Learn More About:** • Does not specifically address B-factors
as predictors of which residues cause docking failures\
• Does not provide a quantitative framework for residue-level
flexibility prediction\
• The gap this review leaves is exactly what our research question
addresses

## 1.2 Kokh, Wade & Wenzel (2011)

**Citation:** Kokh, D. B., Wade, R. C., & Wenzel, W. (2011). Receptor
flexibility in small-molecule docking calculations. WIREs Computational
Molecular Science, 1(2), 298--314. https://doi.org/10.1002/wcms.29

**Summary / Main Point:** A focused review on how receptor
conformational changes affect docking accuracy. The authors argue that
while protein flexibility is critical for understanding receptor--ligand
binding, treating conformational changes remains a major challenge due
to the enormous conformational space that must be sampled and the need
for scoring functions that account for receptor reorganization energy.
Most methods successfully reproduce docking poses of known ligands, but
their benefit for enrichment, affinity prediction, and large-library
screening is less clear. The paper discusses pitfalls and limitations of
flexible-receptor docking strategies.

**What We Can Use:** • Provides the critical perspective: pose
prediction success ≠ virtual screening success\
• Framework for discussing limitations of current flexible-receptor
methods\
• Supports the argument that a simple, fast predictor (B-factor) is
needed to triage which residues matter

**What to Learn More About:** • Does not explore B-factors or any
crystallographic metric for predicting problematic residues\
• Does not address the residue-level question --- only protein-level
flexibility\
• Published in 2011; methods have advanced since, but the fundamental
challenges remain

## 1.3 Cavasotto & Di Maccio (2024)

**Citation:** Cavasotto, C. N., & Di Maccio, E. (2024). Structure and
dynamics in drug discovery. Nature, 1, 1.
https://doi.org/10.1038/s44386-024-00001-2

**Summary / Main Point:** A 2024 Nature review providing a panoramic
view of protein flexibility in drug discovery. Notes that most docking
tools allow full ligand flexibility but keep the protein fixed or with
only limited flexibility. Systematically discusses sources of binding
site flexibility (side-chain rearrangements, loop motions, domain
movements) and their impact on virtual screening and lead optimization.
Emphasizes integrating experimental data (including B-factor analysis)
with computational methods (MD, ensemble docking) to characterize target
flexibility.

**What We Can Use:** • Most up-to-date review (2024) --- establishes
current state of the field\
• Explicitly discusses B-factor analysis as one tool for characterizing
binding site flexibility\
• Published in Nature --- high authority for framing the introduction\
• CRAAP: 24/25

# Section 2: B-factors and Binding Site Flexibility

## 2.1 Halip et al. (2021) --- BFIbs

**Citation:** Halip, L., Avram, S., & Neanu, C. (2021). The B-factor
index for the binding site (BFIbs) to prioritize crystal protein
structures for docking. Structural Chemistry, 32(4), 1693--1699.
https://doi.org/10.1007/s11224-021-01755-5

**Summary / Main Point:** Proposes the Binding site B-factor Index
(BFIbs), defined as the ratio of binding site atom B-factors to
whole-protein atom B-factors. Through automated docking of 26,019
protein--ligand complexes, BFIbs was found to be the only
crystallographic quality metric significantly correlated with docking
pose RMSD (compared to resolution, DPI, R-free, etc.). When BFIbs \< 1,
most docking results achieve RMSD \< 2Å.

**What We Can Use:** • THE most directly relevant paper --- uses
B-factors to predict docking accuracy at the binding site level\
• BFIbs definition and methodology can be directly applied in our
analysis\
• Large-scale statistical evidence (26,019 complexes) supports the
B-factor--docking accuracy link\
• CRAAP: 22/25

**What to Learn More About:** • BFIbs treats the binding site as a whole
--- does not distinguish individual residue contributions\
• Does not explore cross-docking failures or residue-specific B-factor
thresholds\
• Our research extends this to the residue level

## 2.2 Halip et al. (2025) --- LBI

**Citation:** Halip, L., Avram, S., & Neanu, C. (2025). Ligand B-factor
index (LBI): A novel computational index for optimizing protein-ligand
complex docking performance. Structural Chemistry.

**Summary / Main Point:** Extends BFIbs to propose the Ligand B-factor
Index (LBI), defined as the ratio of median binding site B-factors to
median ligand B-factors. On the CASF-2016 benchmark (285 complexes), LBI
shows moderate correlation with experimental binding affinity (Spearman
ρ ≈ 0.48) and outperforms PBI and crystal resolution in pose prediction.
When LBI is in the 0.8--1.2 range, \>75% of scoring functions
successfully predict the native ligand conformation.

**What We Can Use:** • Demonstrates that B-factor ratios (not absolute
values) are more robust predictors\
• Provides a complementary metric to BFIbs for the analysis\
• CRAAP: 22/25

## 2.3 Zhang et al. (2010) --- B-factor & Solvent Accessibility

**Citation:** Zhang, H., Zhang, T., Chen, K., Shen, S., Ruan, J., &
Kurgan, L. (2010). On the relation between residue flexibility and local
solvent accessibility in proteins. Proteins: Structure, Function, and
Bioinformatics, 78(9), 2114--2130. https://doi.org/10.1002/prot.22375

**Summary / Main Point:** Systematically investigates the relationship
between residue flexibility (measured by B-factor) and relative solvent
accessibility (RSA). Finds that local RSA strongly influences B-factor,
and the flexibility-exposure correlation index is strongly related to
amino acid folding stability. This means B-factor is confounded by
solvent exposure --- a critical consideration when interpreting
B-factors at binding sites.

**What We Can Use:** • Provides evidence for B-factor validity as a
flexibility indicator\
• BUT also reveals a key confound: solvent accessibility affects
B-factor\
• Must be discussed in our limitations section\
• CRAAP: 21/25

# Section 3: Impact of Flexibility on Docking Accuracy

## 3.1 Najmanovich et al. (2000) --- Side-chain Flexibility

**Citation:** Najmanovich, R., Kuttner, J., Sobolev, V., & Edelman, M.
(2000). Side-chain flexibility in proteins upon ligand binding.
Proteins: Structure, Function, and Genetics, 39(3), 261--268.
https://doi.org/10.1002/(SICI)1097-0134(20000515)39:3\<261::AID-PROT117\>3.0.CO;2-K

**Summary / Main Point:** Analyzed 980 apo/holo protein pairs from the
PDB. Found that typically only a few binding-site residues change
side-chain conformation upon ligand binding (≤3 residues in \~85% of
cases). Established an amino acid side-chain flexibility scale: Lys \>
Arg, Gln, Met \> Glu, Ile, Leu \> \... \> Phe (Lys bends 25× more
frequently than Phe). CRITICAL FINDING: No correlation between backbone
motion and side-chain flexibility --- meaning B-factors (which primarily
reflect backbone atom motion) may not directly predict side-chain
flexibility.

**What We Can Use:** • Amino acid flexibility scale provides prior
knowledge for predicting problematic residues\
• The \"backbone--side-chain decoupling\" finding is the BIGGEST
challenge to our hypothesis\
• Must be directly addressed: we argue that B-factors capture local
environment dynamics, not just backbone motion per se\
• CRAAP: 21/25

## 3.2 Sherman et al. (2006) --- Induced Fit Docking

**Citation:** Sherman, W., Day, T., Jacobson, M. P., Friesner, R. A., &
Farid, R. (2006). Novel procedure for modeling ligand/receptor induced
fit effects. Journal of Medicinal Chemistry, 49(2), 534--553.
https://doi.org/10.1021/jm050540c

**Summary / Main Point:** Developed the Induced Fit Docking (IFD)
protocol. Tested on 21 drug-relevant cases: rigid docking yielded
average RMSD of 5.5Å (essentially failed), while IFD reduced average
RMSD to 1.4Å (18/21 cases with RMSD ≤ 1.8Å). Provides quantitative proof
that induced fit effects at binding-site residues are the primary cause
of rigid docking failures.

**What We Can Use:** • Quantitative evidence: rigid docking RMSD 5.5Å →
IFD RMSD 1.4Å\
• IFD can serve as a comparison method: does B-factor-guided residue
selection achieve similar improvements?\
• 21-case test set可作为benchmark\
• CRAAP: 23/25

## 3.3 Jain (2009) --- Pocket Adaptation

**Citation:** Jain, A. N. (2009). Effects of protein conformation in
docking: improved pose prediction through protein pocket adaptation.
Journal of Computer-Aided Molecular Design, 23(6), 355--374.
https://doi.org/10.1007/s10822-009-9266-3

**Summary / Main Point:** Quantitatively assessed protein conformation
effects on pose prediction using Surflex-Dock. Cognate docking success
rate: \>75%. Cross-docking with a single fixed protein conformation:
\~20% success. Using multiple conformations with pocket adaptation:
61--75% average success. Demonstrates that binding site conformational
variation is the primary bottleneck for cross-docking.

**What We Can Use:** • Quantitative cross-docking baseline: 20% → 75%
success with multiple conformations\
• Supports the argument that protein conformation is the bottleneck, not
the docking algorithm\
• CRAAP: 24/25

## 3.4 Zhao & Sanner (2007) --- FLIPDock

**Citation:** Zhao, Y., & Sanner, M. F. (2007). FLIPDock: Docking
flexible ligands into flexible receptors. Proteins: Structure, Function,
and Bioinformatics, 68, 726--737. https://doi.org/10.1002/prot.21423

**Summary / Main Point:** Systematically analyzed cross-docking failures
across 4 receptors (cAPK, CDK2, Ricin, HIVp). Of 33 cross-docking
experiments, 25 failed --- and the vast majority of failures were caused
by side-chain conformational changes near the active site. Allowing
these side chains to adopt rotameric conformations successfully rescued
19 complexes (RMSD \< 2.0Å). Using geometric constraints from
bonding/non-bonding interaction networks, 22/25 complexes were
successfully docked.

**What We Can Use:** • Direct evidence: cross-docking failures are
primarily caused by binding-site side-chain flexibility\
• Provides a method for identifying which side chains need to be
flexible (spatial overlap analysis)\
• Our research question: can B-factors replace this more labor-intensive
analysis?\
• CRAAP: 21/25

## 3.5 Gunasekaran & Nussinov (2007) --- Flexible vs Rigid Binding Sites

**Citation:** Gunasekaran, K., & Nussinov, R. (2007). How different are
structurally flexible and rigid binding sites? Sequence and structural
features discriminating proteins that do and do not undergo
conformational change upon ligand binding. Journal of Molecular Biology,
365(1), 257--273. https://doi.org/10.1016/j.jmb.2006.10.038

**Summary / Main Point:** Systematically distinguishes sequence and
structural features of flexible vs rigid binding sites. Flexible binding
sites tend to have higher B-factors, more loop regions, and lower
residue conservation. These features can partially predict whether a
protein will undergo conformational change upon ligand binding.

**What We Can Use:** • B-factor is identified as a KEY distinguishing
feature between flexible and rigid binding sites\
• Provides theoretical support for using B-factors to identify
problematic residues\
• Published in JMB by Nussinov --- high authority\
• CRAAP: 23/25

# Section 4: Ensemble Docking and Flexible-Receptor Methods

## 4.1 Amaro et al. (2018) --- Ensemble Docking Review

**Citation:** Amaro, R. E., Baudry, J., Chodera, J., Demir, Ö.,
McCammon, J. A., Miao, Y., & Smith, J. C. (2018). Ensemble docking in
drug discovery. Biophysical Journal, 114(10), 2271--2278.
https://doi.org/10.1016/j.bpj.2018.02.038

**Summary / Main Point:** Comprehensive review of ensemble docking in
drug discovery. Ensemble docking generates multiple protein
conformations (typically from MD simulations) and docks against each to
approximate receptor flexibility. Discusses conformational sampling
methods (MD, NMR, crystal structure ensembles), conformation selection
strategies (clustering, energy weighting), and applications in virtual
screening.

**What We Can Use:** • Complete methodological framework for ensemble
docking\
• B-factors could serve as a conformation selection criterion ---
discussed but not explored\
• CRAAP: 24/25

## 4.2 Bottegoni et al. (2009) --- 4D Docking

**Citation:** Bottegoni, G., Kufareva, I., Totrov, M., & Abagyan, R.
(2009). Four-dimensional docking: A fast and accurate account of
discrete receptor flexibility in ligand docking. Journal of Medicinal
Chemistry, 52(2), 397--406. https://doi.org/10.1021/jm8009958

**Summary / Main Point:** Developed 4D docking, integrating multiple
receptor conformations into a single docking simulation. Benchmarked on
99 drug-related proteins, 300 ligands, 1113 receptor structures: 77.3%
success rate in correctly reproducing ligand binding geometry. 4× faster
than traditional ensemble docking while maintaining accuracy.

**What We Can Use:** • Large-scale benchmark data for evaluating
B-factor-based methods\
• Demonstrates that crystal structure ensembles (not just MD) can
effectively characterize flexibility\
• CRAAP: 22/25

## 4.3 Carlson (2002) --- Moving Target

**Citation:** Carlson, H. A. (2002). Protein flexibility and drug
design: how to hit a moving target. Current Opinion in Chemical Biology,
6(4), 447--452. https://doi.org/10.1016/S1367-5931(02)00350-4

**Summary / Main Point:** A seminal review that first systematically
articulated the challenge of protein flexibility for drug design.
Introduced the \"moving target\" metaphor --- protein binding sites are
dynamic, not static. Argued that the most advanced drug design and
database mining methods must incorporate protein flexibility.

**What We Can Use:** • Classic conceptual framework for the
introduction\
• \"Moving target\" metaphor is effective for presentations\
• CRAAP: 21/25

## 4.4 De Paris et al. (2018) --- Ensemble Optimizer

**Citation:** De Paris, R., Vahl Quevedo, C., Ruiz, D. D., Gargano, F.,
& de Souza, O. N. (2018). A selective method for optimizing ensemble
docking-based experiments on an InhA fully-flexible receptor model.
BMC Bioinformatics, 19, 235. https://doi.org/10.1186/s12859-018-2228-9

**Summary / Main Point:** Proposes a hierarchical clustering strategy
for ensemble selection using docking-score feedback. Instead of selecting
conformations by structural similarity alone, this method iteratively
refines the ensemble based on how well each conformation performs in
actual docking runs. Tested on InhA (enoyl-ACP reductase) from
*Mycobacterium tuberculosis* with a fully-flexible receptor model.
Demonstrates that score-based ensemble pruning consistently outperforms
RMSD-based clustering for selecting productive receptor conformations.

**What We Can Use:** • Direct methodological competitor: iterative
score-based optimization vs. our B-factor-based selection\
• Serves as a literature benchmark (upper bound) for ensemble
optimization performance\
• Demonstrates that docking-score feedback beats structural similarity
for ensemble selection --- our BFIbs approach is a middle ground\
• CRAAP: 22/25

**What to Learn More About:** • InhA-specific; generalizability to
kinases (CDK2) not tested\
• Requires running docking to select the ensemble --- computationally
more expensive than B-factor pre-filtering\
• Does not leverage crystallographic data (B-factors) at all

## 4.5 Bolstad & Anderson (2009) --- Ensemble Pruning

**Citation:** Bolstad, A. C., & Anderson, A. C. (2009). In pursuit of
virtual lead optimization: Pruning ensembles of receptor structures for
increased efficiency and accuracy during docking. Journal of
Computer-Aided Molecular Design, 23(11), 755--763.
https://doi.org/10.1007/s10822-009-9220-9

**Summary / Main Point:** Introduces a relative-distance-based pruning
method that retains ensemble members whose active-site core geometry
best preserves the conserved ligand-binding pharmacophore. Tested on
dihydrofolate reductase (DHFR): pruned ensembles of 3--5 structures
matched or exceeded the accuracy of full 20+ member ensembles at a
fraction of the computational cost. Key insight: geometric preservation
of the binding site core is more important than conformational diversity
per se.

**What We Can Use:** • Validates the principle that small, well-chosen
ensembles can match large ensembles --- core rationale for BFIbs
selection\
• Geometric pruning provides a complementary benchmark to De Paris
docking-score method\
• Supports our hypothesis that selection strategy matters more than
ensemble size\
• CRAAP: 21/25

**What to Learn More About:** • DHFR-specific; kinase applicability
unknown\
• Geometric pruning uses RMSD of active-site residues --- B-factor may
capture similar information at zero computational cost\
• Published 2009; methods like De Paris (2018) have since advanced the
score-based approach

## 4.6 Bhatt et al. (2024) --- Ensemble Optimizer (EnOpt)

**Citation:** Bhatt, R., Wang, A., & Durrant, J. D. (2024). Teaching
old docks new tricks with machine learning enhanced ensemble docking.
Scientific Reports, 14, 22489.
https://doi.org/10.1038/s41598-024-71699-3

**Summary / Main Point:** Introduces Ensemble Optimizer (EnOpt), a
machine-learning tool that improves the accuracy and interpretability of
ensemble virtual screening. EnOpt trains interpretable ML models on
per-ligand docking scores across multiple receptor conformations to
predict which conformations are most likely to produce active ligands.
Unlike prior ensemble selection methods (De Paris 2018, Bolstad 2009),
EnOpt learns from the docking data itself rather than relying on
structural clustering or geometric pruning. Demonstrated improved
enrichment over traditional averaging and heuristic selection across
multiple targets.

**What We Can Use:** • STATE-OF-THE-ART ensemble optimizer --- this is
the primary literature benchmark for our ensemble selection comparison\
• ML-based approach provides the upper bound against which BFIbs-guided
selection should be compared\
• Interpretable ML framework: feature importance reveals which
conformations drive performance --- analogous to our B-factor analysis\
• Directly addresses the same problem we tackle: which structures to
include in an ensemble\
• CRAAP: 23/25

**What to Learn More About:** • Requires running full docking before
selection (like De Paris) --- BFIbs is pre-docking and zero-cost\
• ML training overhead may be impractical for small-scale studies\
• Does not leverage crystallographic B-factors --- our approach is
complementary, not competitive\

## 4.7 Bhatt et al. (2023) --- EnOpt (Biophysical Journal)

**Citation:** Bhatt, R., Wang, A., & Durrant, J. D. (2023). Ensemble
optimizer: Interpretable scoring functions for virtual screening.
Biophysical Journal, 122(3), 67a.
https://doi.org/10.1016/j.bpj.2022.11.557

**Summary / Main Point:** Earlier conference presentation of the
Ensemble Optimizer concept. Proposes interpretable scoring functions for
combining and comparing ensemble docking scores, moving beyond simple
averaging heuristics. Establishes the mathematical framework later
expanded in the 2024 Scientific Reports paper.

**What We Can Use:** • Provides the theoretical foundation for EnOpt\
• Demonstrates that heuristic averaging (common in ensemble docking) is
suboptimal --- validates our motivation for principled selection\
• CRAAP: 20/25 (conference abstract, less detail than full paper)

# Summary: Literature by Priority

## Must-Read (Core Support)

  ----------------------------------------------------------------------------
  **Paper**         **Score**         **Key Contribution**   **Role in
                                                             Argument**
  ----------------- ----------------- ---------------------- -----------------
  Halip et al.      22/25             B-factor → docking     Core evidence for
  (2021) BFIbs                        accuracy (26K          hypothesis
                                      complexes)             

  Sherman et al.    23/25             Rigid 5.5Å → IFD 1.4Å  Quantifies the
  (2006) IFD                          RMSD                   problem

  Najmanovich et    21/25             Backbone--side-chain   Key challenge to
  al. (2000)                          decoupling             address

  Zhao & Sanner     21/25             25/33 failures from    Direct evidence
  (2007) FLIPDock                     side-chain changes     

  Jain (2009)       24/25             Cross-docking 20% →    Quantitative
                                      75%                    baseline

  Gunasekaran &     23/25             B-factor distinguishes Theoretical
  Nussinov (2007)                     flexible/rigid sites   support

  Lexa & Carlson    N/A               Complete methods       Methodological
  (2012)                              landscape              framework
  ----------------------------------------------------------------------------

## Important References

  -----------------------------------------------------------------------
  **Paper**         **Score**         **Key             **Role in
                                      Contribution**    Argument**
  ----------------- ----------------- ----------------- -----------------
  Cavasotto & Di    24/25             Latest Nature     Current state of
  Maccio (2024)                       review on         field
                                      flexibility in    
                                      drug discovery    

  Kokh et al.       N/A               Pose ≠ screening; Critical
  (2011)                              limitations of    perspective
                                      flexible docking  

  Halip et al.      22/25             B-factor ratio    Complementary
  (2025) LBI                          metric for pose   method
                                      prediction        

  Zhang et al.      21/25             B-factor          Limitation to
  (2010)                              confounded by     discuss
                                      solvent           
                                      accessibility     

  Amaro et al.      24/25             Ensemble docking  Context for
  (2018)                              framework         solutions

  Bottegoni et al.  22/25             4D docking        Benchmark data
  (2009)                              benchmark (99     
                                      proteins)         

  De Paris et al.   22/25             Score-based       Literature
  (2018)                              ensemble          benchmark for
                                      optimizer;         ensemble
                                      iterative         optimization
                                      docking feedback  

  Bolstad &         21/25             Geometric         Validates
  Anderson (2009)                     pruning retains   small-ensemble
                                      pharmacophore;    principle
                                      3--5 ≈ 20+        

  Bhatt et al.      23/25             ML-based          State-of-the-art
  (2024) EnOpt                        Ensemble          benchmark for
                                      Optimizer;        ensemble
                                      interpretable     selection
                                      scoring           

  Carlson (2002)    21/25             \"Moving target\" Introduction
                                      conceptual        framing
                                      framework         
  -----------------------------------------------------------------------

# Identified Research Gap

The literature reveals a clear gap that our research question addresses:

- What is known:

<!-- -->

- Binding site flexibility is the primary cause of cross-docking
  failures (Jain 2009; FLIPDock; Sherman 2006)

- B-factors correlate with overall binding site docking quality (BFIbs,
  Halip 2021)

- Flexible and rigid binding sites can be distinguished by features
  including B-factor (Gunasekaran & Nussinov 2007)

<!-- -->

- What remains unsolved:

<!-- -->

- No study has validated at the RESIDUE level whether B-factors predict
  which specific flexible residues cause docking failures

- Najmanovich (2000) found backbone motion and side-chain flexibility
  are uncorrelated --- this challenges using B-factors to predict
  side-chain behavior, but has not been systematically tested in the
  docking context

- No study has compared B-factor-guided residue selection with existing
  methods (IFD, ensemble docking) in terms of efficiency and accuracy
