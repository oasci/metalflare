#import "lib.typ": config, todo, addcite

///////////////
// Variables //
///////////////

#let TITLE = "Unraveling the Cu(I) Sensing Mechanism of roGFP2"
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
    authors: (),
    version: VERSION,
    bibliography: bibliography(BIBTEX_PATH, title: "References", style: "american-chemical-society")
)

= Methods

== Molecular simulations

Four simulation sets---reduced, oxidized, Cu(I), and Na#super[+]---were performed to gain atomistic insight into roGFP2 sensing mechanisms.
The reduced simulations represent our designed control system with protonated Cys147 and Cys204.
We will refer to this crucial cysteine pair as the 'sensing cysteines'.
Oxidized simulations contain a disulfide bond linking the sensing cysteines.
Cu(I) simulations let the ion coordinate with the sensing cysteines and use their respective metal-coordinated parameters.
Na#super[+] simulations use the same Cys metal-coordinated parameters without adding additional ions, which offers a control for Cys parameters and ion identity.

=== Protein preparation

Initial protein structures for reduced (#link("https://files.rcsb.org/download/1JC0.pdb")[1JC0]) and oxidized (#link("https://files.rcsb.org/download/1JC1.pdb")[1JC1]) states of roGFP2 @hanson2004investigating were retrieved from the Protein Data Bank (PDB). @berman2000protein
The structures were processed using in-house Python and bash scripts (available free of charge at #link("github.com/oasci/metalflare")[github.com/oasci/metalflare]).
The first chain, along with the crystallographic water molecules, were centered to the origin and rotated to minimize the box size using NumPy,@harris2020array SciPy,@virtanen2020scipy and MDAnalysis @michaud2011mdanalysis packages.

Given the relevance and availability of the GFP mechanism's force field parameters, the chromophore was modeled in its anionic state.
All selenomethionine residues (MSE) were converted to methionine (MET), and Cys147 and Cys204 were transformed into the appropriate states.
Glu222 was protonated to model the I state of the GFP photocycle; this still probes the anionic chromophore's stability while offering insight into ground-state proton transfer (GSPT) dynamics.
The protonation states of all other residues were determined with PDB2PQR @jurrus2018improvements, using the default parameters.
The pdb4amber tool was then used to validate the PDB file before proceeding.

=== Simulation preparation

System preparation was performed using the tleap module of AmberTools v23.6.
The protein structure was parameterized using the ff19SB force field. @tian2019ff19sb
For the solvent environment, we employed the OPC3 water model, @izadi2016accuracy known for its balanced representation of water properties in biomolecular simulations.
The 12-6 nonbonded model and parameters for all ions were taken from Sengupta et al. @sengupta2021parameterization
The system was neutralized by adding Na#super[+] and Cl#super[-] ions as needed.
Additional ions were introduced to achieve a solvent ionic strength of 0.150 M to mimic physiological conditions.
The protein was solvated in a rectangular box, with a minimum distance of 10 $angstrom$ between the protein and the box edges to minimize periodic boundary condition artifacts.
Parameters from Breyfogle et al. @breyfogle2023molecular were employed for the anionic chromophore.

=== Minimization

The prepared system underwent a four-stage energy minimization protocol using Amber23 #addcite() to relieve any unfavorable interactions and optimize the structure.
All minimization stages used the steepest descent method for the first 1000 steps, followed by the conjugate gradient method for the remaining steps, with a maximum of 5000 steps per stage.
A non-bonded cutoff of 10.0 $angstrom$ was applied throughout.
Periodic boundary conditions were employed, and coordinates were wrapped to the primary unit cell.
The minimization progress was monitored by writing energies every step and coordinates every 200 steps.

#emph[Stage 1:] Initial minimization was performed with restraints (force constant: 5.0 kcal/mol/$angstrom$#super[2]) on all non-hydrogen atoms of the entire system, allowing hydrogen atoms to relax and adjust their positions.
#emph[Stage 2:] The system was further minimized with restraints (same force constant) on all non-hydrogen atoms of the solute, excluding water molecules and ions, allowing solvent and ions to equilibrate around the solute.
#emph[Stage 3:] Minimization continued with reduced restraints (force constant: 2.0 kcal/mol/$angstrom$#super[2]) applied only to the protein backbone, allowing side chains and other flexible parts to relax.
#emph[Stage 4:] Final minimization was performed with further reduced restraints (force constant: 1.0 kcal/mol/$angstrom$#super[2]).
The resulting minimized structure was the starting point for subsequent relaxation and production simulations.

=== Relaxation

Following energy minimization, the system underwent a three-stage relaxation protocol using Amber23 to gradually equilibrate the structure and solvent.
Three independent runs were initiated with random initial velocities drawn from the Maxwell-Boltzmann distribution at 100 K.
All subsequent simulations were continued using the respective run's restart files.

#emph[Stage 1:] An initial 20 ps NVT (constant Number of particles, Volume, and Temperature) simulation was performed with a 2 fs time step. The system was heated from 100 K to 300 K using Langevin dynamics with a collision frequency of 5 ps⁻¹. Restraints (force constant: 1.0 kcal/mol/$angstrom$#super[2]) were applied to the protein backbone.
SHAKE algorithm was used to constrain bonds involving hydrogen atoms.
The non-bonded cutoff was set to 10.0 $angstrom$.
#emph[Stage 2:] A 1 ns NPT simulation followed, maintaining the temperature at 300 K using Langevin dynamics (collision frequency: 5 ps⁻¹).
The pressure was regulated at 1.01325 bar using the Monte Carlo barostat with a relaxation time of 1 ps.
Restraints on the same atoms were reduced (force constant: 0.5 kcal/mol/$angstrom$#super[2]).
#emph[Stage 3:] The final relaxation stage consisted of a 1 ns NPT simulation with all positional restraints removed.

=== Production simulations

All production runs were performed under the same setup as the last relaxation stage.
Each run was simulated for 500 ns with coordinates saved every 10 ps.
The resulting trajectories from all three replicates were used for subsequent analyses, providing a substantial 1.5 $mu$s of simulation data for each system state.

=== Hydrogen bonding cutoff

A hydrogen bond of X---H $dots.c$ Y---Z, where X is the donor, and Y is the acceptor atom, can be classified based on distances and angles.
One characteristic recommended by IUPAC is that the H $dots.c$ Y distance is less than the sum of H and Y van der Waals radii.
Hydrogen (1.10 $angstrom$) and oxygen (1.52 $angstrom$) @mantina2009consistent would have a cutoff of 2.62 $angstrom$.
Others @hubbard2010hydrogen recommend a cutoff of 2.50 $angstrom$ based on structural analysis @mcdonald1994satisfying and quantum chemical calculations @liu2008geometrical.
Since the difference between a 2.5 and 2.62 $angstrom$ cutoff is likely a substantially weak hydrogen bond, we will use an H $dots.c$ Y cutoff of 2.5 $angstrom$.

= Results

Our molecular dynamics simulations of roGFP2 in its reduced, oxidized, and Cu(I)-bound states reveal significant conformational changes.
Oxidized and reduced simulations closely match experimental structures (PDB IDs: #link("https://www.rcsb.org/structure/1JC0")[1JC0] and #link("https://www.rcsb.org/structure/1JC1")[1JC1]), with mean C#sub[$alpha$]-C#sub[$alpha$] distance differences less than 0.04 $angstrom$ (@fig-md\B).
An experimental structure of roGFP2-Cu(I) is currently unavailable; however, our simulations suggest that it is 0.44 $angstrom$ larger than the reduced state.
Cu(I) coordinated between Cys147 and Cys204 throughout all simulations, enhancing conformational flexibility.
For example, @fig-md\B shows a broader distribution of C#sub[$alpha$]-C#sub[$alpha$] distances with a slight bimodal peak centered around 4.48 $angstrom$ and 4.96 $angstrom$.

Notably, the hydrogen bonding behavior of Thr203 shown in @fig-md\C exhibits striking alterations across the different states.
Thr203 shows dynamic interactions with the chromophore in the reduced state, evidenced by two primary distance peaks at 1.75 and 5.32 $angstrom$.
Oxidation disrupts this hydrogen bond, shifting the dominant peak to 5.83 $angstrom$ and decreasing the hydrogen bonding population.
Conversely, Cu(I) binding enhances the Thr203--chromophore hydrogen bond, resulting in a sharp peak at 1.75 $angstrom$---more than tripling the corresponding population in the reduced state.
These findings suggest that Cu(I) induces a conformational change that optimally positions Thr203 for chromophore interaction.

#figure(
    image(FIG_DIR + "shared-figures/001-md-overview/fig001.png", width: 6.0in),
    caption: [
        Analysis of molecular dynamics simulations of roGFP2 in four Cys147 and 204 states: reduced, oxidized, Cu(I)-bound, and Na#super[+]-bound.
        *(A):* A visual depiction of key residues surrounding the chromophore in roGFP2.
        Cys147 and Cys204 are not depicted but are adjacent to His148 and Thr203, respectively.
        *(B):* Probability density function for the C#sub[$alpha$]-C#sub[$alpha$] distance between Cys147 and Cys204.
        *(C):* Distance between Thr203 hydroxyl (HG1) and Cro66 oxygen (OH) as a function of state, showing how Cu(I) induces substantially more Thr203 stabilization of the anionic chromophore.
        *(D):* Distance between Glu222 (HE2) and Ser205 (OG), representing changes in the ground-state proton transfer pathway.
        Cu(I) binding and oxidation disrupt the interaction critical for chromophore reprotonation.
    ],
    placement: auto
) <fig-md>

Further analysis of the ground-state proton transfer (GSPT) pathway reveals that oxidation or Cu(I) binding could specifically inhibit the proton transfer from Glu222 to Ser205 (@fig-md\D).
In the reduced state, this initial step occurs with a probability of 0.299, supporting a functional GSPT pathway.
However, the probability drops significantly in the oxidized (0.000) and Cu(I)-bound (1.122$times$10#super[-4]) states, as Glu222 hydrogen bonds with the chromophore instead of Ser205.
Despite this inhibition, the probabilities for proton transfer from Ser205 to water and from water to the chromophore remain comparable across all states.
These results highlight that oxidation and Cu(I) binding could induce specific conformational changes that disrupt the initial GSPT step.

= Discussion

Our molecular simulations suggest that the oxidation of Cys147 and 204 primarily destabilizes the anionic chromophore by eliminating Thr203's hydrogen bonding.
While His148 and Try145 coordinate more with the chromophore, the dynamics are vastly similar to reduced and Cu(I) states, as shown in @fig-pes-his-tyr.
These conditions may favor the stabilization of the neutral chromophore, leading to an increasing A-band, as experimentally observed.
The GSPT pathway between Glu222 and Ser205 was disrupted in the oxidized and Cu(I) states.
This could impede the reprotonation of the anionic chromophore, but it is not a differentiating factor between the two states.

In the case of Cu(I) binding, our simulations indicate destabilization of key structural elements---such as the $beta$-sheet between His148 and Thr203---and reduced hydrogen bonding through His148 and Tyr145.
These changes may facilitate non-radiative decay pathways, selectively quenching the B-band fluorescence without significantly affecting the A-band.
The strong coordination of Thr203 with the chromophore's phenolate oxygen in the Cu(I)-bound state might help maintain the anionic chromophore population, unlike the oxidized state, yet the overall structural alterations promote quenching.

Na#super[+] simulations aimed to discern whether the effects observed are unique to Cu(I) or merely a consequence of nonspecific cation binding.
The results revealed that the Cys147--Cys204 C#sub[$alpha$] distance in the Na#super[+] simulations was similar to that of the Cu(I)-bound state, although the peak was slightly longer at 5.06 $angstrom$.
Crucially, the Na#super[+] simulations did not exhibit enhanced hydrogen bonding between Thr203 and the chromophore observed in the Cu(I)-bound state.
Instead, the Thr203--chromophore interactions in the Na#super[+] simulations resembled those of the reduced state.
Additionally, the GSPT pathway dynamics in the presence of Na#super[+] closely mirrored those of the reduced state, with a functional initial proton transfer step from Glu222 to Ser205.
These findings suggest that the structural alterations and functional disruptions seen with Cu(I) binding are not replicated with Na#super[+].
This underscores the specificity of Cu(I) in inducing conformational changes that lead to selective quenching of the anionic chromophore fluorescence.

We must interpret these findings cautiously due to the limitations, including the lack of classical force field parameters for the excited and neutral chromophore states.
Our data still offers atomistic insight that correlates with the experimental observations.
Further studies---incorporating excited-state dynamics and explicit modeling of the neutral chromophore---are imperative to unravel the complex sensing mechanisms of roGFP2 fully.

#pagebreak()
= Supplemental information

== Cys147 and Cys207 C#sub[$alpha$] distance

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

== $beta$-strand fraying

GFP and its variants are characterized by a distinctive $beta$-barrel structure, which plays a crucial role in their fluorescence properties.
The structural integrity of this $beta$-barrel is essential for maintaining the protein's fluorescence characteristics and its sensitivity to environmental changes.
In roGFP2, Cys147 is strategically positioned near the C-terminus of a $beta$-pleated sheet, flanked by His148 and Thr203.

#figure(
    image(FIG_DIR + "f-backbone/f001-his148_h-thr203_o/f001-his148_h-thr203_o-pdf.svg", width: 3.5in),
    caption: [
        Probability density of His148 and Thr203 backbone hydrogen bonding under various roGFP2 conditions.
    ],
    placement: bottom
) <fig-beta-his147>
Our molecular dynamics simulations reveal a striking structural change upon Cu(I) binding to roGFP2.
Specifically, we observe a significant disruption in the $beta$-sheet structure near Cys147.
@fig-beta-his147 illustrates that Cu(I) binding leads to the breaking of a key $beta$-sheet hydrogen bond in this region.
This observation is further substantiated by quantitative analysis presented in @tab-beta-his147, which shows a dramatic decrease in the hydrogen bond probability between the backbone --NH and C=O groups within a 2.5 $angstrom$ distance.
The probability drops from 0.865 in the reduced state to a mere 0.063 in the Cu(I)-bound state, indicating a near-complete loss of this stabilizing interaction.

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

The structural perturbation at this site could also propagate through the protein structure.
Indeed, the intermolecular distance peak between Asn146 and Ser205 increases from 3.80 to 4.48 $angstrom$ with bound Cu(I).
Distances with Cu(I) extend to at least 0.66 $angstrom$ further than the reduced state.
Under oxidizing conditions, the peak distance only decreases to 3.56 $angstrom$ and hydrogen bonds 14.4 % during the simulations.

No substantial changes were observed in the Val150 and Leu201 $beta$-strand.

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

== Disruption of ground state proton transfer

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

