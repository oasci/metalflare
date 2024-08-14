#import "lib.typ": config, todo, addcite

///////////////
// Variables //
///////////////

#let TITLE = "Unraveling the roGFP2 Cu(I) Sensing Mechanism"
#let VERSION = (
  "v" + datetime.today().display("[year].[month padding:none].[day padding:none]")
)
#let BASE_DIR = "../../../"
#let BIBTEX_PATH = BASE_DIR + "management/literature/refs.bib"
#let FIG_DIR = BASE_DIR + "figures/"

/////////////
// Content //
/////////////

#show: config.with(
    title: [#TITLE],
    authors: (
        (
        name: "Alex M. Maldonado",
        department: [Department of Biological Sciences],
        organization: [University of Pittsburgh],
        location: [Pennsylvania, United States],
        email: "alex.maldonado@pitt.edu"
        ),
    ),
    version: VERSION,
    bibliography: bibliography(BIBTEX_PATH, title: "References", style: "american-chemical-society")
)

Here, we propose an atomistic mechanism for Cu(I) sensing for a redox-sensitive green fluorescent protein variant, roGFP2.
#todo("This is a very rough draft and not polished. Ideas and concepts are stable.")

= GFP fluorescence mechanism

The fluorescence mechanism in green fluorescent protein (GFP) is a complex interplay of photophysical and photochemical processes that occur at the molecular level.
This section delves into the current understanding of fluorescence in GFP, exploring the fundamental principles of excitation and de-excitation, as well as the various factors that influence these processes.
We will examine how the protein environment modulates the fluorescent properties, the critical role of chromophore protonation states, and the intricate dynamics of excited-state phenomena such as non-adiabatic crossings and proton transfer.
From this point forward, we will refer to enhanced GFP (eGFP) as "GFP" and eGFP chromophore @cormack1996facs as "chromophore" that is the result of the F64L and S65T mutations from the #emph[Aequorea victoria] wild type GFP (wtGFP) @prasher1992primary.

== Chromophore

#figure(
    image(FIG_DIR + "h-background/h005-cro/cro-a.svg", width: 2.3in),
    caption: [
        Neutral (i.e., A state) GFP chromophore.
    ],
    placement: auto
) <fig-cro-a>
Fluorescence in all GFPs begins with the absorption of a photon by the chromophore, typically in the blue region of the visible spectrum (\~300 to 500 nm).
This absorption process is linked to the unique molecular structure of the chromophore.
The chromophore from eGFP is formed autocatalytically from three amino acid residues---Thr65, Tyr66, and Gly67---through a series of cyclization, dehydration, and oxidation reactions.
The resulting structure consists of a hydroxybenzylidene imidazolinone moiety, which forms an extended $pi$-conjugated system, as shown in @fig-cro-a.
wtGFP chromophore has a primary alcohol group (i.e., Ser65) instead of a secondary (i.e., Thr65).

=== Excitation

The $pi$-conjugated system is a cornerstone of the chromophore's light-absorbing properties.
Photons of specific energies, precisely corresponding to the energy gap between the ground state ($S_0$) and the chromophore's first excited state ($S_1$), excite the delocalized electrons in the conjugated bonds.
Planarity, protonation state, and, notably, the protein environment collectively contribute to fine-tuning the exact absorption wavelength.

When a blue photon is absorbed, it promotes an electron from the highest occupied molecular orbital (HOMO) to the lowest unoccupied molecular orbital (LUMO) of the chromophore.
This electronic transition is predominantly $pi → pi^*$, reflecting the excitation within the $pi$-conjugated system.

=== Protonation states

The chromophore protonation state is a critical determinant of its photophysical properties.
The chromophore can exist in two primary forms: a neutral (protonated) and anionic (deprotonated) state, each exhibiting distinct spectroscopic signatures.

#figure(
    image(FIG_DIR + "h-background/h005-cro/cro-b.svg", width: 2.3in),
    caption: [
        Anionic (i.e., B-state) GFP chromophore.
    ],
    placement: auto
) <fig-cro-b>
In its neutral form (shown in @fig-cro-a), the chromophore's phenolic oxygen is protonated, resulting in an absorption maximum typically around 395-400 nm.
This is often referred to as the 'A state'.
The anionic form of the chromophore, where the phenolic oxygen is deprotonated, is primarily responsible for GFP's characteristic green fluorescence.
This state, often called the B-state, has an absorption maximum of approximately 475-490 nm.

Several factors influence the relative population of these two states, including the local pH, specific interactions within the protein environment, and mutations in the protein sequence.

=== Planarity

Chromophore planarity plays a pivotal role in determining its fluorescent properties, particularly through enhanced conjugation and increased quantum yield.
Enhanced conjugation in a planar chromophore structure is primarily due to the maximized overlap of p-orbitals in the $pi$-conjugated system.
This optimal alignment of p-orbitals, perpendicular to the molecular plane, facilitates efficient delocalization of $pi$-electrons across the entire conjugated system.
Such extensive electron delocalization has profound effects on the chromophore's electronic structure, most notably in reducing the HOMO-LUMO gap.

Planarity significantly impacts the chromophore's quantum yield by restricting non-radiative decay pathways.
In non-planar configurations, rotation around single bonds, can serve as an efficient route for non-radiative decay.
However, when the chromophore is held in a planar conformation by the protein environment, these rotational movements are severely limited.
This restriction of molecular motion creates energy barriers in the excited state potential energy surface, effectively preventing the chromophore from accessing geometries that favor non-radiative decay.
Consequently, with fewer available non-radiative pathways, the excited chromophore is more likely to return to the ground state via fluorescence emission, directly increasing the fluorescence quantum yield.

=== Neutral-chromophore fluorescence

Excitation and subsequent emission in the A state without structural changes is one (infrequent) possibility.

#figure(
    image(FIG_DIR + "h-background/h006-cro-excitation/gfp-a-emission.png", width: 2.0in),
    caption: [
        Excitation and emission of neutral chromophore.
    ],
    placement: auto
) <fig-a-emission>
Excitation of the neutral state occurs at approximately 395 nm, corresponding to the absorption of violet-blue light.
This excitation promotes the chromophore to $S_1$ without immediate proton transfer.
The subsequent emission from this excited neutral state results in weak blue fluorescence with a peak around 460 nm.

This protonated chromophore can be further stabilized by hydrogen bonding with a water molecule, Thr203, or Ser205.
The presence or absence of these interactions contribute to the distinct excitation and emission.
Although, this emission is relatively weak compared to the anionic state discussed next.
#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-a-other.png", width: 3.0in),
    caption: [
        Example stabalizing configuration for the neutral chromophore.
    ],
    placement: auto
) <fig-a-residues>

=== B-band fluorescence

The anionic state of the chromophore represents a key configuration responsible for the protein's characteristic green fluorescence.
In this state, the chromophore exists in its deprotonated form, with the phenolic oxygen carrying a negative charge.

Excitation of the anionic chromophore occurs at approximately 475 nm, corresponding to the absorption of blue light.
This excitation promotes the chromophore to its first excited singlet state ($S_1$).
The subsequent relaxation and emission result in the bright green fluorescence typically associated with GFP, with a peak around 508 nm.

#figure(
    image(FIG_DIR + "h-background/h006-cro-excitation/gfp-b-emission.png", width: 2.0in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
) <fig-b-emission>

Specific interactions within the protein barrel stabilize the anionic state.
Notably, Thr203 plays a crucial role in stabilizing the anionic form through hydrogen bonding with the deprotonated phenolic oxygen.
A protonated Glu222 could also form a hydrogen bond with the anionic chromophore, further contributing to its stabilization.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-b-other-nolabel.png", width: 3.0in),
    caption: [
        Example stabalizing configuration for the anionic chromophore.
    ],
    placement: auto
) <fig-b2-emission>

The protein environment around the chromophore is critical in maintaining this anionic configuration.
The $beta$-barrel structure of GFP provides a hydrophobic pocket that shields the chromophore from bulk solvent, contributing to the high quantum yield of fluorescence in this state.

=== Excited-state proton transfer

The excited-state proton transfer (ESPT) mechanism represents a fundamental process in GFP, contributing significantly to its unique spectroscopic properties.
This process involves the initial excitation of the neutral chromophore, followed by rapid, successive proton transfer events in the excited state, ultimately resulting in emission from an anionic species.
An overview of the process is shown below and is referred to as 'A band'.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-photocycle.png", width: 3.5in),
    caption: [
        The GFP photocycle involved an excited- and ground-state proton transfer.
    ],
    placement: auto
) <fig-espt-emission>

1. Upon absorption of a photon at approximately 395 nm, the neutral chromophore is promoted to $S_1$.

2. In this excited state, the chromophore exhibits markedly different acid-base properties compared to its ground state, becoming a much stronger acid.
    This enhanced acidity facilitates the transfer of a proton from the chromophore to a proximal acceptor within the protein matrix.
    The proton transfer pathway involves a sophisticated hydrogen-bonding network.
    A critical component of this network is a strategically positioned water molecule, which acts as the initial proton acceptor.
    This water molecule is part of a proton wire that includes Ser205 and terminates at Glu222.
    The transfer occurs on an ultrafast timescale, typically in the order of picoseconds.

3. Following the ESPT, the system exists transiently in an intermediate state (I\*), characterized by an anionic chromophore and a protonated Glu222.
    I\* subsequently relaxes and emits fluorescence at approximately 508 nm, closely resembling the B-state emission profile.

4. The final step is a reverse ground-state proton transfer (GSPT) from Glu222 through Ser205, a water molecule, and terminated at the chromophore.

Note that while the chromophore is protonated when excited, this fluorescence is a separate pathway from the 'neutral-chromophore fluorescence' mechanism discussed above.

== Roles of other crucial residues

The exceptional fluorescent properties of GFP arise from a complex interplay of molecular interactions within its $beta$-barrel structure.

=== His148

Structurally, His148 helps maintain $beta$-strand stability by hydrogen bonding to Arg168's backbone #addcite().
His148 directly interacts with the chromophore, forming a critical hydrogen bond.
Specifically, the imidazole side chain of His148 donates a hydrogen bond to the chromophore's phenolate oxygen #addcite().
This interaction aids in locking the chromophore into a planar conformation and minimizes non-radiative decay pathways.
Upon chromophore excitation, His148 can also help facilitate proton movement @shinobu2010visualizing.
Others have also explored introducing a hole in the $beta$ barrel through an H148G mutation, allowing metal ions to diffuse into the protein and interact directly with the chromophore @barondeau2002structural.

=== Thr203

Thr203 is located close to the chromophore within the $beta$-barrel structure.
Its hydroxyl group can hydrogen bond with the anionic chromophore's phenolate oxygen in specific conformations.
Mutating Thr203 to Ile or Val destabilizes the B state and dramatically reduces its excitation peak, thus driving the ground-state equilibrium towards the A state @kummer2000effects.
Yellow FP (YFP) is also derived from mutations by placing chemical groups that can $pi$--$pi$-stacking with the chromophore @zimmer2002green.
Quenching of A\* emissions through an ESPT is still present in T203 mutants; however, T203I and T203V mutants are slower down while T203Y maintains wild-type speeds @jung2005photophysics.

=== Tyr145

Tyrosine's hydroxyl group produces a bulge in the $beta$-barrel, and a phenylalanine mutation enhances the protein's thermal stability @nasu2021structure.
This mutation alleviates local structural strain, thereby increasing the thermal resilience of the GFP without disrupting its overall $beta$-barrel fold @akiyama2012experimental.
To the best of our knowledge, Tyr145 hydrogen bonding to the chromophore does not play a large role in B-state stabilization---especially since phenylalanie does not significantly affect the B band.

= roGFP2 contains redox-sensing cysteines

Redox-sensitive GFPs (roGFPs) are engineered variants of GFP designed to report cellular redox states through changes in their fluorescence properties @hanson2004investigating.
Among these, roGFP2 has emerged as a particularly useful probe due to its ratiometric readout and midpoint potential.
It is suitable for measuring redox conditions in reducing cellular compartments like the cytosol and mitochondria.

roGFP2 contains two cysteine residues (S147C and Q204C mutations) on adjacent $beta$-strands near the chromophore of eGFP.
(A structural depiction of the relevant residues is in @fig-rogfp2-structure.)
The strategically placed cysteines can form a reversible disulfide bond in response to changes in the surrounding redox environment.
When a disulfide bond is formed, it alters the preferred chromophore protonation state in roGFP2.
This alteration leads to reciprocal changes in the excitation peaks, which shift to approximately 400 nm and 490 nm.

#figure(
    image(FIG_DIR + "h-background/h009-rogfp2-sims/relevant-reduced.png", width: 3.0in),
    caption: [
        Fluorescence-relevant residues in reduced roGFP2.
    ],
    placement: auto
) <fig-rogfp2-structure>

= roGFP2 fluorescence response

@fig-rogfp2-fluorescence illustrates the redox-dependent fluorescence properties of roGFP2, showing the protein's excitation spectra under various redox conditions.
These spectra were obtained by monitoring emission at 511 nm while scanning excitation wavelengths from 350 to 500 nm.
The protein samples were equilibrated in buffers containing different ratios of oxidized and reduced dithiothreitol (DTT) to achieve a range of defined redox potentials.

roGFP2 exhibits two distinct excitation peaks: (1) the A band at 400 nm and (2) the B band at 490 nm.
These peaks correspond to 511 nm fluorescence of the chromophore's excited neutral (protonated) and anionic (deprotonated) forms.
We will refer to the neutral and anionic forms of the chromophore as "A-state" and "B-state", respectively.

#figure(
    image(FIG_DIR + "h-background/h002-rogfp2-fluorescence/rogfp2-redox-fluorescence.png", width: 3.0in),
    caption: [
        Relative fluorescence at 511 nm after excitation scan from 350 to 500 nm at -0.310 (blue), -0.275 (green), and -0.240 V (orange) redox potentials.
        Adapted with permission from Hanson et al. @hanson2004investigating distributed under the CC-BY-4.0 license.
    ],
    placement: auto
) <fig-rogfp2-fluorescence>

We observe a clear shift in the excitation spectrum as the environment becomes more oxidizing (moving from -0.240 V to -0.310 V).
The intensity of the A band increases significantly, while the B band experiences a concomitant decrease.
This spectral shift reflects the formation of the disulfide bond between Cys147 and Cys204 and an increasing preference for the protonated chromophore.

Notably, the spectra display a clear isosbestic point at approximately 425 nm, where the fluorescence intensity remains constant regardless of the redox state.
This isosbestic point is a hallmark of a two-state system, confirming that roGFP2 is transitioning cleanly between its oxidized and reduced forms without significant intermediate states.
This spectral feature is a critical indicator of the probe's behavior: it suggests that the engineered disulfide bond is forming and breaking as intended, without competing side reactions or alternative conformations significantly affecting the fluorescence.
The maintenance of this isosbestic point across various redox potentials implies that the structural changes induced by oxidation and reduction are consistent and reversible.
Any deviation from this behavior, such as a shift in the isosbestic point or its disappearance, would suggest more complex interactions---potentially involving intermediate states, protein structural changes, or interactions with other molecules---that could complicate data interpretation.

== Possible perturbations

Changes in A- or B-band emissions could indicate various environmental changes in the chromophore.

- *Equilibrium ratio of A- or B-state populations.*
    A- or B-state chromophore stability changes impact the relative proportion of A- and B-band absorbance.
    For example, roGFP1 has a larger A band than B band @hanson2004investigating.
    A-state stability is directly correlated to A-band fluorescence, B-band stability would correlate to B-band fluorescence, etc.
- *Excited-state proton transfer (ESPT) from A\* $arrow.r$ I\*.*
    Emissions at the typical 511 nm (green) fluorescence from the A-state require an ESPT from Cro66 to Glu222 through a coordinated water molecule and Ser205.
    Prohibiting ESPT would result in more radiative emission at ~460 nm.
- *Ground-state proton transfer (GSPT) from I $arrow.r$ A.*
    Reprotonating the chromophore through a GSPT is crucial for maintaining the A- and B-band lifetime.
    Disrupting the Glu222 $arrow.r$ Ser205, Ser205 $arrow.r$ H#sub[2]O, or H#sub[2]O $arrow.r$ Cro66 pathway would decrease the A-state population---likely with a corresponding B-state increase.
- *Non-radiative emissions.*
    Enhanced flexibility of the chromophore through protein conformational shifts would lead to additional non-radiative emissions, thereby lowering the relative fluorescence in that state.
- *Changes in chromophore environment.*
    Alterations in the local environment around the chromophore can significantly impact its fluorescence characteristics. These changes may arise from various cellular processes or experimental conditions.
    For example, shifts in the polarity of the chromophore's immediate surroundings can affect the stability of its protonated (A) and deprotonated (B) states. Increased exposure to water molecules might preferentially stabilize the deprotonated form, potentially enhancing B-band fluorescence relative to the A-band.
    Molecular species capable of quenching fluorescence can dramatically affect roGFP2's signal intensity if they come close to the chromophore. Common cellular quenchers include molecular oxygen and specific metal ions like Cu#super[2+] or Fe#super[3+].
- *Protein structural changes.*
    While the formation of the disulfide bond between Cys147 and Cys204 is the primary structural change driving roGFP2's redox response, other conformational changes in the protein could also impact its fluorescence properties:
    Increased flexibility of the $beta$-barrel structure enclosing the chromophore might allow for more non-radiative decay pathways. This could reduce overall fluorescence intensity or changes in the relative intensities of the A and B bands.


= Molecular simulations

Four simulation sets, reduced, oxidized, Cu(I), and Na#super[+], were performed to gain atomistic insight into roGFP2 sensing mechanisms.
The reduced simulations represent our designed control system with protonated Cys147 and Cys204.
We will refer to this crucial cysteine pair as the 'sensing cysteines', a key focus of our research.
Oxidized simulations contain a disulfide bond linking the sensing cysteines.
Cu(I) simulations let the ion coordinate with the sensing cysteines and use their respective metal-coordinated parameters.
Na#super[+] simulations use the same metal-coordinated parameters of Cu(I) simulations without adding additional ions, which offers a control for the ion identity and cysteine parameters.

== Protein preparation

Initial protein structures for reduced (#link("https://files.rcsb.org/download/1JC0.pdb")[1JC0]) and oxidized (#link("https://files.rcsb.org/download/1JC1.pdb")[1JC1]) states of roGFP2 were retrieved from the Protein Data Bank (PDB).
The structures were processed using in-house Python and bash scripts (available free of charge at #link("https://gitlab.com/oasci/studies/metalflare")[gitlab.com/oasci/studies/metalflare]).
The first chain, along with the crystallographic water molecules, were centered to the origin while minimize the box values using NumPy #addcite(), SciPy #addcite(), and MDAnalysis #addcite() packages.

Given the relevance and availability of the GFP mechanism's force field parameters, the chromophore was modeled in its anionic state.
All selenomethionine residues (MSE) were converted to methionine (MET), and Cys147 and Cys204 were transformed into the appropriate residues.
Glu222 was protonated to model the I state of the GFP photocycle; this still probes the anionic chromophore's stability while offering insight into GSPT dynamics.
The protonation states of all other residues were determined with PDB2PQR @jurrus2018improvements, using the default parameters to ensure standardization.
The pdb4amber tool was then used to validate the PDB file before proceeding.

== Cu(I) placement

#todo("Add GFN2-xTB calculations.")

== Simulation preparation

System preparation was performed using the tleap module of AmberTools v23.6 #addcite().
The protein structure was parameterized using the ff19SB force field #addcite().
For the solvent environment, we employed the OPC3 water model, #addcite(), known for its balanced representation of water properties in biomolecular simulations.
The 12-6 nonbonded model and parameters for all ions were taken from Sengupta et al. @sengupta2021parameterization
The system was neutralized by adding Na+ and Cl- ions as needed.
Additional ions were introduced to achieve a solvent ionic strength of 0.150 M to mimic physiological conditions.
The protein was solvated in a rectangular box, with a minimum distance of 10 $angstrom$ between the protein and the box edges to minimize periodic boundary condition artifacts.
Parameters from Breyfogle et al. @breyfogle2023molecular were employed for the anionic chromophore.

== Minimization

The prepared system underwent a four-stage energy minimization protocol using Amber23 #addcite() to relieve any unfavorable interactions and optimize the structure.
All minimization stages used the steepest descent method for the first 1000 steps, followed by the conjugate gradient method for the remaining steps, with a maximum of 5000 steps per stage.
A non-bonded cutoff of 10.0 $angstrom$ was applied throughout.
Periodic boundary conditions were employed, and coordinates were wrapped to the primary unit cell.
The minimization progress was monitored by writing energies every step and coordinates every 200 steps.

#emph[Stage 1:] Initial minimization was performed with restraints (force constant: 5.0 kcal/mol/$angstrom$#super[2]) on all non-hydrogen atoms of the entire system, allowing hydrogen atoms to relax and adjust their positions.
#emph[Stage 2:] The system was further minimized with restraints (same force constant) on all non-hydrogen atoms of the solute excluding water molecules and ions, allowing solvent and ions to equilibrate around the solute.
#emph[Stage 3:] Minimization continued with reduced restraints (force constant: 2.0 kcal/mol/$angstrom$#super[2]) applied only to the protein backbone, allowing side chains and other flexible parts to relax.
#emph[Stage 4:] Final minimization was performed with further reduced restraints (force constant: 1.0 kcal/mol/$angstrom$#super[2]).
The resulting minimized structure served as the starting point for subsequent relaxation and production simulations.

== Relaxation simulations

Following energy minimization, the system underwent a three-stage relaxation protocol using Amber23 to equilibrate the structure and solvent gradually.
Three independent runs were initiated with random initial velocities to ensure adequate sampling.
All subsequent simulations were continued using the respective run's restart files.

#emph[Stage 1:] An initial 20 ps NVT (constant Number of particles, Volume, and Temperature) simulation was performed with a 2 fs time step. The system was heated from 100 K to 300 K using Langevin dynamics with a collision frequency of 5 ps⁻¹. Restraints (force constant: 1.0 kcal/mol/$angstrom$#super[2]) were applied to the protein backbone.
SHAKE algorithm was used to constrain bonds involving hydrogen atoms.
The non-bonded cutoff was set to 10.0 $angstrom$.
#emph[Stage 2:] A 1 ns NPT simulation followed, maintaining the temperature at 300 K using Langevin dynamics (collision frequency: 5 ps⁻¹).
Pressure was regulated at 1.01325 bar using the Monte Carlo barostat with a relaxation time of 1 ps.
Restraints on the same atoms were reduced (force constant: 0.5 kcal/mol/$angstrom$#super[2]).
#emph[Stage 3:] The final relaxation stage consisted of a 1 ns NPT simulation with all positional restraints removed.

== Production simulations

All production runs were performed under the same setup as the last relaxation stage.
Each run was simulated for 500 ns with coordinates saved every ten ps.
The resulting trajectories from all three replicates were used for subsequent analyses, providing a cumulative 1.5 $mu$s of simulation data for the system under study.


== Analysis

From the MD trajectories, we extracted two primary sets of data:

1. Dihedral angles ($phi$, $psi$) were calculated using MDAnalysis #addcite() for specific residues known to be crucial for roGFP2 function.
2. Key interactions between various residues through intermolecular distances.

Dihedral angles were transformed using

$ frac(1 - cos (theta), 2) $

This transformation maps the circular dihedral data to a [0, 1] range, preserving the periodicity while differentiating between cis and trans conformations.

All input features ($X$) were standardized using sklearn's `StandardScaler` to ensure each feature contributes equally to the model.
This transformation ensures continuity across the periodic boundary, facilitating more accurate modeling.

=== Coordinated water detection

#todo("Write this.")

=== Potential of mean force (PMF)

#todo("Write this.")

=== Feature correlation

To investigate the relationship between structural descriptors and various features, we employed Partial Least Squares (PLS) regression analysis.
This multivariate statistical technique was chosen for its ability to handle high-dimensional, correlated data and reveal underlying patterns in complex datasets.

A PLS regression model was fitted to $X$ and response variable $y$ using sklearn's `PLSRegression` with two components for each simulation state.
The model's performance was evaluated using the $R^2$ score.

Data points were projected onto the space of the first two PLS components.
A 2D histogram was created in this space, with bin colors representing the mean $y$.

Loading vectors for each feature were plotted as arrows in the PLS component space.
The magnitude and direction of these arrows indicate the importance and relationship of each feature to the PLS components.
The magnitude of each feature's loading vector was calculated as the Euclidean norm of its first two PLS components.

A dashed line representing the direction of maximum change in the response variable was added to the plot.
This line, referred to as the derivative line, indicates the direction in the PLS component space along which $y$ increases most rapidly.

== Feature importance

We developed a machine learning pipeline to elucidate the relationship between backbone dihedral angles and the target feature.
Our approach employed two complementary models:

- *XGBoost Regressor:* A gradient boosting algorithm chosen for its high performance and ability to capture non-linear relationships #addcite().
- *Elastic Net:* A regularized linear regression model selected to identify linear correlations while mitigating multicollinearity #addcite().

The dataset was randomly partitioned into training (80%) and testing (20%) sets, ensuring model generalizability.
We utilized `GridSearchCV` with 3-fold cross-validation to tune model hyperparameters.
For XGBoost, we optimized the number of estimators (250-700), learning rate (0.05-0.2), maximum tree depth (5-9), and regularization terms ($alpha$: 0.0-0.2, $gamma$: 0.8-1.0).
For Elastic Net, we tuned the regularization strength ($alpha$: 10#super[-5] - 5) and L1 ratio (0.2-1.0).
Performance was assessed using mean squared error (MSE) and coefficient of determination ($R^2$) on the held-out test set.

== Hydrogen bond cutoff

A hydrogen bond of X---H $dots.c$ Y---Z, where X is the donor and Y is the acceptor atom, can be classified based on distances and angles.
One characteristic recommended by IUPAC is that the H $dots.c$ Y distance is less than the sum of H and Y van der Waals radii.
Hydrogen (1.10 $angstrom$) and oxygen (1.52 $angstrom$) @mantina2009consistent would have a cutoff of 2.62 $angstrom$.
Others @hubbard2010hydrogen recommend a cutoff of 2.50 $angstrom$ based on structural analysis @mcdonald1994satisfying and quantum chemical calculations @liu2008geometrical.
Since the difference between a 2.5 and 2.62 $angstrom$ cutoff is likely a substantially weak hydrogen bond, we will use a H $dots.c$ Y cutoff of 2.5 $angstrom$.

= Fluorescence mechanism of Cu(I) distinct from oxidation

Joel Rosenbaum performed in vitro assays of roGFP2 under H#sub[2]O#sub[2] (i.e., oxidation) and Cu(I) conditions by monitoring 528 nm emissions #todo("Check emission value") with scanning excitation wavelengths between 380 and 500 nm.
Relative fluorescence with respect to #todo("How was this normalized? Check methods") is shown in @fig-cu-fluorescence.

#figure(
    image("s1b.png", width: 3.5in),
    caption: [
        Relative fluorescence of roGFP2 under reduced, oxidized, and Cu(I) conditions from 380 to 500 nm excitation scan.
        #todo("Check emission value?")
        Apo (i.e., reduced) roGFP2 exhibits typical bimodal absorption of A-band (excited at 400 nm) and B-band (excited at 488) peaks.
        Upon roGFP2 oxidation from 1 mM H#sub[2]O#sub[2], a shift in A- (increased) and B-band (decreased) absorption and subsequent 528 nm emission marks a corresponding change in neutral and anionic chromophore populations.
        Binding of Cu(I), however, exhibits a larger decrease in B-band without the A-band increase observed when oxidized.
    ],
    placement: auto
) <fig-cu-fluorescence>

Apo represents our control with roGFP2 in its reduced state.
Adding 1.0 mM H#sub[2]O#sub[2] results in the expected A- and B-band proportions shift.
Again, this indicates that oxidizing Cys147 and Cys204 and forming a disulfide bridge increases the neutral chromophore equilibrium proportion---with a corresponding decrease in the anionic state.
This oxidation-induced change in the chromophore protonation state is the basis for roGFP2's utility as a redox sensor in biological systems.

Cu(I), on the other hand, exhibits an entirely different fluorescence response.
Strikingly, a mere 1 $mu$M of Cu(I) drastically reduces the B-band intensity while keeping the A band relatively stable.
The distinct responses to oxidation and Cu(I) binding provide compelling evidence that Cu(I) binding to roGFP2 does not simply alter the chromophore's protonation state preference, as oxidation does. Instead, these results suggest that Cu(I) may disturb the anionic chromophore's electronic state or induce conformational changes that specifically affect its environment, underscoring the need for further exploration.

= Cu(I) binding enhances roGFP2 backbone flexibility

First, we investigate the backbone dynamics over 1.5 $mu$s of MD simulations.

== Cys147 and Cys207 C#sub[$alpha$] distance

Experimental structures of both the reduced (PDB ID: #link("https://www.rcsb.org/structure/1JC0")[1JC0]) and oxidized (PDB ID: #link("https://www.rcsb.org/structure/1JC1")[1JC1]) states of for roGFP2 exhibited a mean C#sub[$alpha$]-C#sub[$alpha$] distance of 4.30 ± 0.12 and 4.07 ± 0.09 $angstrom$, respectively @hanson2004investigating.
Our MD simulations agreed well with experimental observations as shown in @tab-alpha-c with mean differences less than 0.04 $angstrom$.

#figure(
    table(
        columns: (auto, auto, auto),
        stroke: none,
        [*State*], [*Experimental ($angstrom$)*], [*Simulations ($angstrom$)*],
        table.hline(),
        [*Reduced*], [4.30 ± 0.12], [4.34 ± 0.47],
        [*Oxidized*], [4.07 ± 0.09], [4.11 ± 0.29],
        [*Cu(I)*], [N/A], [4.78 ± 0.82],
        [*Na#super[+]*], [N/A], [5.06 ± 0.89],
    ),
    placement: none,
    caption: [Mean and standard deviations of C#sub[$alpha$]-C#sub[$alpha$] distances between Cys147 and Cys204.]
) <tab-alpha-c>

The standard deviations, $sigma$, also provide an indication of conformational flexibility.
Unsurprisingly, the $sigma$ of oxidized roGFP2 exhibits is substantially smaller than the reduced state.

An experimental structure of roGFP2-Cu(I) is currently not available; however, our simulations suggest that it is on average 0.44 $angstrom$ larger than the reduced state.
Cu(I) stayed coordinated between Cys147 and Cys204 throughout all simulations while enhancing conformational flexibility.
For example, @fig-alpha-c shows a broader distribution of C#sub[$alpha$]-C#sub[$alpha$] distances with a slight bimodal peaks centered around 4.48 $angstrom$ and 4.96 $angstrom$.

#figure(
    image(FIG_DIR + "b-cys/b004-cys147_ca-cys204_ca/b004-cys147_ca-cys204_ca-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
) <fig-alpha-c>

== $beta$-strand fraying

GFP and its variants are characterized by a distinctive $beta$-barrel structure, which plays a crucial role in their fluorescence properties.
This $beta$-barrel is composed of 11 $beta$-strands, forming a robust scaffolding that encapsulates and protects the centrally located chromophore.
The structural integrity of this $beta$-barrel is essential for maintaining the protein's fluorescence characteristics and its sensitivity to environmental changes.
In roGFP2, Cys147 is strategically positioned near the C-terminus of a $beta$-pleated sheet, flanked by His148 and Thr203 (Figure @fig-rogfp2-structure).

#figure(
    image(FIG_DIR + "f-backbone/f001-his148_h-thr203_o/f001-his148_h-thr203_o-pdf.svg", width: 3.5in),
    caption: [
        Probability density of His148 and Thr203 backbone hydrogen bonding under various roGFP2 conditions.
    ],
    placement: auto
) <fig-beta-his147>
Our molecular dynamics simulations reveal a striking structural change upon Cu(I) binding to roGFP2.
Specifically, we observe a significant disruption in the $beta$-sheet structure near Cys147.
@fig-beta-his147 illustrates that Cu(I) binding leads to the breaking of a key $beta$-sheet hydrogen bond in this region.
This observation is further substantiated by quantitative analysis presented in @tab-beta-his147, which shows a dramatic decrease in the hydrogen bond probability between the backbone --NH and C=O groups within a 2.5 $angstrom$ distance.
The probability drops from 0.865 in the reduced state to a mere 0.063 in the Cu(I)-bound state, indicating a near-complete loss of this stabilizing interaction.
(Unsurprisingly, oxidized roGFP2 stabilizes this hydrogen bond.)

This disruption of the $beta$-sheet hydrogen bonding network can be characterized as "fraying" of the $beta$-strand.
Fraying is a phenomenon where the regular hydrogen bonding pattern at the termini of secondary structure elements, particularly $beta$-strands, becomes disrupted, leading to increased local flexibility.

#figure(
    caption: [Hydrogen bonding probability between residue backbones],
    table(
        columns: (auto, auto, auto),
        stroke: none,
        [*State*], [*His148 - Thr203*], [*Asn146 - Ser205*],
        table.hline(),
        [*Reduced*], [0.865], [0.047],
        [*Oxidized*], [0.997], [0.144],
        [*Cu(I)*], [0.063], [0.000],
        [*Na#super[+]*], [0.380], [0.002],
    ),
    placement: auto
) <tab-beta-his147>

In the case of roGFP2, the Cu(I)-induced fraying at the C-terminus of this $beta$-strand is likely to have several important implications.
Given the proximity of this structural change to the chromophore, it may directly influence the chromophore's electronic environment, contributing to the observed changes in fluorescence properties upon Cu(I) binding.

#figure(
    image(FIG_DIR + "f-backbone/f002-ser205_h-asn146_o/f002-ser205_h-asn146_o-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
) <fig-beta-his146>
The structural perturbation at this site could also propagate through the protein structure.
Indeed, @fig-beta-his146 shows the intermolecular distance peak between Asn146 and Ser205 increasing from 3.80 to 4.48 $angstrom$ with bound Cu(I).
Distances with Cu(I) extend to at least 0.66 $angstrom$ further than the reduced state.
Under oxidizing conditions, the peak distance only decreases to 3.56 $angstrom$ and hydrogen bonds 14.4 % during the simulations.

No changes are observed in the Val150 and Leu201 $beta$-strand.

= Cu(I) binding affects anionic chromophore hydrogen binding

As previously mentioned, several residues interact with the chromophore.
Numerous studies have observed complex interactions that have led to several variants tailored for specific applications.
roGFP2 is no different.

== His148 and Tyr145

Our analysis of hydrogen bonding interactions with the chromophore in roGFP2 revealed distinct patterns across different protein states.
@fig-pes-his-tyr illustrates the 2D potential energy surfaces for Tyr145 and His148 interactions with the chromophore.
#figure(
    image(FIG_DIR + "g-cro-interact/g001-his148_hd1-tyr145_hh/g001-pes-combined.png", width: 3.5in),
    caption: [
        2D potential of mean force (PMF) plot of the distance between the chromophore's phenolate oxygen and His148 (x-axis) and Tyr145 (y-axis).
        The color map applies a linear colormap normalized from 0 (blue) to 4 (yellow) kcal/mol.
    ],
    placement: auto
) <fig-pes-his-tyr>

In the reduced state of roGFP2, we observed a diverse energy landscape characterized by multiple local minima.
A minimum is located at Tyr145-chromophore and His148-chromophore distances of approximately 1.8 $angstrom$ and 1.85 $angstrom$, respectively.
This minimum corresponds to a configuration where both residues simultaneously form hydrogen bonds with the chromophore.
Two other local minima were observed, representing configurations where Tyr145 or His148 independently formed a hydrogen bond with the chromophore.
Still, both residues exhibit conformational flexibility where neither residues stabilize the anionic chromophore.
This pattern suggests significant conformational flexibility in the reduced state, allowing for various hydrogen bonding arrangements.
#figure(
    caption: [Structural hydrogen bonding probability to chromophore],
    table(
        columns: (auto, auto, auto, auto, auto,),
        stroke: none,
        [*State*], [*Thr203*], [*His148*], [*Tyr145*], [*Glu222*],
        table.hline(),
        [*Reduced*], [0.191], [0.486], [0.612], [0.689],
        [*Oxidized*], [0.009], [0.691], [0.795], [0.977],
        [*Cu(I)*], [0.619], [0.339], [0.641], [0.895],
        [*Na#super[+]*], [0.299], [0.347], [0.671], [0.627],
    ),
    placement: auto
) <tab-cro-hbond>

Upon oxidation of roGFP2, we noted a marked change in the energy landscape.
The oxidized state exhibited a single, pronounced global minimum at similar Tyr145-chromophore (1.80 $angstrom$) and His148-chromophore (1.91 $angstrom$) distances as the reduced state.
However, this energy well was deeper and narrower than the reduced state.
The absence of significant additional local minima indicates a strong preference for the configuration where both Tyr145 and His148 simultaneously hydrogen bond with the chromophore.
This suggests that oxidation induces a more rigid and stable arrangement of these critical residues around the chromophore.

By contrast, the Cu(I)-bound state of roGFP2 exhibited similar minima as reduced simulations but with a notable increase in flexibility with shallower minima.
The dual-hydrogen bonding minimum was significantly less pronounced than the reduced and oxidized states.
It is worth noting that Tyr145 is 1.9 times more likely to hydrogen bond to the chromophore than His148, indicating an asymmetry in their hydrogen bonding behavior.

His148 and Thr203 $beta$-strand fraying appears to correlates correlate with decreased stabilization of the anionic chromophore through His148.
#todo("Add origin analysis if needed")

== Thr203

Our molecular dynamics simulations reveal a striking change in the hydrogen bonding behavior of Thr203.
In the reduced state of roGFP2, we observe several distribution peaks of Thr203-chromophore distances.
Two primary peaks at 5.32 and 1.75 $angstrom$ suggest dynamic populations across the simulations.
#figure(
    image(FIG_DIR + "g-cro-interact/g010-cro66_oh-thr203_hg1/g010-cro66_oh-thr203_hg1-pdf.svg", width: 3.5in),
    caption: [
        Probability density of the distance between Thr203 HG and the chromophore's phenolate oxygen.
        The grey region indicates our hydrogen-bonding region.
    ],
    placement: auto
)

Oxidation of roGFP2 dramatically alters this interaction.
The hydrogen bonding population decreases substantially, with a minor peak at 1.80 $angstrom$.
Instead, we observe a dominant peak at 5.83 $angstrom$, indicating that oxidation largely disrupts the Thr203-chromophore hydrogen bond.
This disruption likely contributes to the shift in chromophore protonation state observed upon oxidation.

Remarkably, Cu(I) binding to roGFP2 significantly enhances the Thr203-chromophore hydrogen bond.
We observe a sharp, dominant peak at 1.75 $angstrom$, which is more than three times the corresponding peak in the reduced state.
This dramatic increase in hydrogen bonding suggests that Cu(I) binding induces a conformational change that positions Thr203 optimally for interaction with the chromophore.

#todo("Add origin analysis if needed")

= Disruption of ground-state proton transfer

Our simulations reveal significant changes in the GSPT pathway of roGFP2 under different conditions.
The GSPT pathway---which involves proton transfer from Glu222 to Ser205, a water molecule, and finally the chromophore---shows distinct behavior in reduced, oxidized, and Cu(I)-bound states of roGFP2.
@tab-gspt-hbond presents the observed probabilities of the intermolecular distance to be less than 2.5 $angstrom$.
#figure(
    caption: [GSPT step probability],
    table(
        columns: (auto, auto, auto, auto),
        stroke: none,
        [*State*], [*Glu222 $arrow.r$ Ser205*], [*Ser205 $arrow.r$ H#sub[2]O*], [*H#sub[2]O $arrow.r$ Cro66*],
        table.hline(),
        [*Reduced*], [0.299], [0.516], [0.560],
        [*Oxidized*], [0.000], [0.416], [0.601],
        [*Cu(I)*], [1.122$times$10#super[-4]], [0.517], [0.538],
        [*Na#super[+]*], [0.295], [0.472], [0.404],
    ),
    placement: auto
) <tab-gspt-hbond>

In the reduced state of roGFP2, we observe a probability of 0.299 for the initial proton transfer step from Glu222 to Ser205.
The subsequent steps in the pathway, from Ser205 to water and from water to Cro66, occur with probabilities of 0.516 and 0.465, respectively.
These values suggest a functional GSPT pathway in the reduced state, with each step occurring at a substantial frequency throughout the simulations.

Strikingly, we observe a complete inhibition of the initial GSPT step in oxidized simulations and a rare occurrence in Cu(I); instead, Glu222 hydrogen bonds to the secondary alcohol group of the chromophore instead of Ser205.
This represents a marked change from the reduced state and indicates a significant alteration in the protein's proton transfer capabilities upon oxidation or Cu(I) binding.

Despite the inhibition of the initial step, the probabilities for the latter two steps of the GSPT pathway remain relatively consistent across all three states.
In the oxidized state, the probabilities of the Ser205 to water and water to Cro66 steps are 0.416 and 0.379, respectively.
Similarly, in the Cu(I)-bound state, these steps show probabilities of 0.517 and 0.452.
These values are comparable to those observed in the reduced state, suggesting that the capability for these latter proton transfer steps is maintained even when the initial step is inhibited.

#figure(
    image(FIG_DIR + "e-proton-wire/e001-ser205_og-glu222_he2/e001-ser205_og-glu222_he2-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
 ],
    placement: auto
) <fig-glu-ser>
Further examination of the Glu222-Ser205 interaction reveals intriguing details about the GSPT pathway dynamics.
Figure @fig-glu-ser illustrates the probability density of the distance between Glu222 and Ser205 across the three states of roGFP1.

We observe a bimodal distribution of the Glu222-Ser205 distance in the reduced state.
The primary peak at approximately 1.81 $angstrom$ corresponds to a hydrogen bonding distance.
A secondary peak at around 4.54 $angstrom$ suggests an alternative conformation in which Glu222 does not interact directly with Ser205.
Interestingly, this bimodal distribution in the reduced state indicates that Glu222 can alternate between interacting with Ser205 and the chromophore.

In stark contrast, the oxidized and Cu(I)-bound states show an unimodal distribution centered around 3.93 and 4.74 $angstrom$, respectively.
This single peak at a larger distance indicates that Glu222 and Ser205 are consistently separated in these states and do not form a direct hydrogen bond.
This observation directly corresponds to the zero probability of proton transfer from Glu222 to Ser205 in both oxidized and Cu(I)-bound states.

These findings collectively demonstrate that oxidation and Cu(I) binding induce specific conformational changes in roGFP2 that increase the Glu222-Ser205 distance, thereby inhibiting the initial step of the GSPT pathway.
This inhibition occurs without significantly altering the latter steps of the pathway, pointing to a localized effect on the Glu222-Ser205 interaction rather than a global disruption of the proton transfer network.

= Proposed Cu(I) sensing mechanism

In the fluorescence experiments, we observe the following characteristics.

- *Oxidized*:
 - A-band increases.
 - B-band decreases.
 - An isosbestic point indicates there is an equilibrium between the two states.
- *Cu(I)*:
 - A-band remains relatively stable but will begin to decrease with increasing Cu(I) concentration.
 - B-band dramatically decreases even at low Cu(I) concentrations.
 - Lacks an isosbestic point.
 - Neutral and anionic chromophores still absorb at expected wavelengths (data not shown).

These results---agreeing with previous studies #addcite()---indicate that oxidizing Cys147 and Cys204 leaves the fluorescence mechanism intact.
However, the neutral chromophore becomes slightly preferred.

Cu(I)-binding appears to quench the B-band primarily.
The lack of the isosbestic point and no A-band increase indicates no detectable shift in chromophore protonation preference.
However, there is likely a population shift along with fortuitous A-band quenching that makes it appear unchanged.
Applying Occam's razor would discourage this mechanism and instead attribute Cu(I) sensing to quenching the B-band faster than the A-band.

Our molecular simulations suggest the following concerning the reduced state.

- *Oxidized*:
 - Added strain between Cys147 and Cys204.
 - Stabilized $beta$-sheet between His148 and Thr203.
 - Enhanced chromophore stabilization through His148 and Tyr145 hydrogen bonding.
 - Thr203 rarely coordinates with the chromophore.
 - Potentially broken GSPT pathway to re-protonate the chromophore.
- *Cu(I)*:
 - Additional backbone flexibility around Cys147 and Cys204.
 - Destabilizes $beta$-sheet between His148 and Thr203.
 - Less hydrogen bonding of anionic chromophore through His148 and Tyr145.
 - Thr203 strongly favors coordinating with the phenolate oxygen in the chromophore.
 - Potentially broken GSPT pathway to re-protonate the chromophore.

Based on the experimental evidence and simulation data, we propose a multi-faceted atomistic mechanism for Cu(I) sensing in roGFP2.
In the reduced state, roGFP2 maintains a dynamic equilibrium between the neutral and anionic chromophore states, with the anionic state being favored.
Thr203 plays a crucial role in stabilizing the anionic chromophore through hydrogen bonding with the phenolate oxygen. The functional GSPT pathway allows for efficient re-protonation of the chromophore, maintaining the A-band fluorescence.

Upon Cu(I) binding, roGFP2 undergoes distinct conformational changes that affect the anionic chromophore state.
The simulations reveal increased backbone flexibility around the Cu(I) binding site, destabilization of the $beta$-sheet structure, and reduced hydrogen bonding of the anionic chromophore through His148 and Tyr145.
These changes likely contribute to the quenching of B-band fluorescence by promoting non-radiative decay pathways.
Despite the destabilizing effects of the conformational changes, the strong coordination of Thr203 to the chromophore's phenolate oxygen in the Cu(I)-bound state may help maintain the anionic chromophore population.
This stabilizing effect of Thr203 could explain the absence of a significant increase in the A-band fluorescence upon Cu(I) binding, as observed in the experimental data.

The dramatic decrease in B-band fluorescence upon Cu(I) binding can be attributed to the quenching of the anionic chromophore through non-radiative decay pathways, which are likely promoted by the destabilization of the $beta$-sheet structure and reduced hydrogen bonding through His148 and Tyr145.
The potentially disrupted GSPT pathway may further contribute to this quenching by impeding the re-protonation of the anionic chromophore.
The relatively stable A-band in the presence of Cu(I) suggests that the GSPT pathway is not completely inhibited, allowing for some re-protonation of the chromophore to maintain the neutral state population.

In summary, we propose a multi-faceted Cu(I) sensing mechanism in roGFP2 that involves several factors.
Thr203 maintains the anionic chromophore population through hydrogen bonding with the phenolate oxygen. The destabilization of the $beta$-sheet structure and reduced hydrogen bonding through His148 and Tyr145 promote non-radiative decay pathways, leading to the quenching of B-band fluorescence.
The partially functional GSPT pathway allows for some re-protonation of the chromophore, maintaining the neutral state population and contributing to the relatively stable A-band fluorescence.

= Limitations

This section details my critical review of this work.
Feedback will help strengthen, adjust, or clarify points in this report.

- *Absence of an Isosbestic Point:*
 The lack of an isosbestic point in the Cu(I) fluorescence data is interpreted as evidence for a different mechanism from oxidation. This interpretation is reasonable, but additional controls should support it.
 Have other metal ions been tested to ensure this response is specific to Cu(I)?
 This would help to exclude the possibility that the observed changes are due to general metal-induced conformational shifts rather than specific Cu(I) interactions.
- *Role of Thr203:*
 The report suggests that Thr203 plays a crucial role in maintaining the anionic chromophore population through hydrogen bonding in the Cu(I)-bound state.
 However, the simulation results indicate a strong preference for Thr203 to coordinate with the chromophore in the presence of Cu(I), contrasting with its behavior in the oxidized state.
 This dual role needs to be reconciled—does Thr203 stabilize the anionic chromophore, or does it contribute to its quenching by altering the electronic environment?
 This point needs clarification, possibly with additional experimental validation (e.g., mutagenesis studies).
- *Disrupted GSPT Pathway:*
 The report proposes that Cu(I) binding disrupts the GSPT pathway, leading to quenching of B-band fluorescence.
 However, the evidence for this disruption is primarily computational.
 Experimental data, such as proton exchange rates or hydrogen bonding network perturbations (e.g., through NMR or FTIR spectroscopy), would strengthen this argument.
 Additionally, the complete inhibition of the GSPT pathway in the oxidized state, as suggested by the simulations, needs to be better supported by the existing literature on GFPs and their variants.
- *Simulation of Cu(I) Binding:*
 The choice of force fields and Cu(I) coordination accuracy need further justification.
 Cu(I) has complex coordination chemistry, and improper parametrization could lead to erroneous results.
 Have other methods, such as quantum mechanical/molecular mechanical (QM/MM) simulations, been considered to validate the Cu(I) coordination observed in the MD simulations?
 This would address potential inaccuracies in the force field's Cu(I) treatment.
- *Backbone Flexibility and $beta$-Strand Fraying:*
 The report links Cu(I) binding to increased backbone flexibility and $beta$-strand fraying, leading to decreased stabilization of the anionic chromophore.
 While the data support this, it remains speculative without direct experimental evidence.
 Could experimental techniques, such as hydrogen-deuterium exchange (HDX) mass spectrometry, corroborate these findings?
- *Hydrogen Bonding Analysis:*
 The hydrogen bonding probabilities and distances presented are critical to the proposed mechanism.
 However, the methodology for calculating these probabilities needs more detail.
 Were these calculated over the entire simulation time, or were specific time windows chosen?
 How were fluctuations in hydrogen bond lengths and angles accounted for in the analysis?
 The stability and significance of these hydrogen bonds need to be better contextualized, perhaps by comparing them to known benchmarks in similar systems.
- *Fluorescence Spectra Interpretation:*
 The interpretation of the fluorescence spectra, especially the distinct responses of the A- and B-bands, relies heavily on the computational predictions.
 It would be beneficial to see more direct experimental tests, such as site-directed mutagenesis of Thr203, His148, or Tyr145, to confirm their roles in Cu(I) sensing.
 How do mutations in these residues affect the fluorescence response to Cu(I)? This would directly test the proposed mechanism.
- *Comparison with Other Metal Ions:*
 To further substantiate the specificity of the Cu(I) sensing mechanism, comparisons with other metal ions (e.g., Zn(II), Fe(III), Ni(II)) should be included.
 Does roGFP2 exhibit a similar fluorescence response to these metals, or is the response unique to Cu(I)?
 This experiment would help to rule out non-specific metal binding as a confounding factor.
- *Mechanistic Plausibility:*
 The proposed Cu(I) sensing mechanism involves multiple factors, such as the disruption of $β$-strand integrity, altered hydrogen bonding, and the inhibition of the GSPT pathway. While this multi-faceted approach is thorough, the sheer complexity of the mechanism raises questions about its biological plausibility. Are all these events necessary and sufficient to explain the observed fluorescence changes? The report could benefit from a discussion on the relative importance of each factor, potentially simplifying the proposed mechanism to focus on the most critical elements. Simplifying the narrative where possible would make the mechanism more digestible and potentially more convincing.
- *Energetic Considerations:*
 The report suggests that Cu(I) binding leads to increased backbone flexibility and destabilization of the $β$-sheet, affecting chromophore stability and fluorescence. However, the energetic cost of these structural rearrangements is not discussed. Given that a delicate balance of forces typically stabilizes protein structures, substantial conformational changes could have significant energetic penalties. It would be useful to discuss whether these conformational changes are likely to occur spontaneously in the presence of Cu(I) or whether additional factors (e.g., cellular environment, protein-protein interactions) might influence this process.
- *Sampling Adequacy:*
 The report mentions a total of 1.5 μs of MD simulation data, which is substantial. However, it is not clear whether this sampling is sufficient to capture the full range of conformational changes and dynamic behavior of roGFP2 in the presence of Cu(I). Have convergence tests been performed to ensure that the observed changes are not artifacts of insufficient sampling? Additionally, are there significant differences between the three independent simulation runs, or are the results consistent across all replicates? A more detailed analysis of the sampling adequacy and reproducibility of the results would strengthen the computational aspect of the study.
- *Force Field Limitations:*
 While the report uses well-established force fields (e.g., ff19SB for proteins), the limitations of these force fields should be acknowledged, especially in the context of metal coordination chemistry. Force fields are often parametrized for standard amino acid residues and may not accurately capture the unique coordination environment of Cu(I). This could lead to artifacts or misinterpretations of the Cu(I) binding effects. A brief discussion on the limitations of the chosen force fields and any potential impact on the results would be beneficial.
- *Normalization of Fluorescence Data:*
 The report mentions that the fluorescence data were normalized, but the normalization method is unclear. Given that the interpretation of the fluorescence spectra is central to the proposed mechanism, the normalization method should be described in detail. Was the fluorescence intensity normalized to the protein concentration, emission wavelength, or another factor? Inconsistent or inappropriate normalization could lead to misleading conclusions, so this aspect needs careful attention.

