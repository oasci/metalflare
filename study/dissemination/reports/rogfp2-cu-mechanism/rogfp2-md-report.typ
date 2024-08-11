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

Here, we propose atomistic mechanisms for a redox-sensitive green fluorescent protein variant roGFP2 including a novel Cu(I) sensing application.

= GFP fluorescence mechanism

The fluorescence mechanism in green fluorescent protein (GFP) is a complex interplay of photophysical and photochemical processes that occur at the molecular level.
This section delves into the current understanding of fluorescence in GFP, exploring the fundamental principles of excitation and de-excitation, as well as the various factors that influence these processes.
We will examine how the protein environment modulates the fluorescent properties, the critical role of chromophore protonation states, and the intricate dynamics of excited-state phenomena such as non-adiabatic crossings and proton transfer.
From this point forward, we will refer to eGFP as "GFP" and eGFP chromophore @cormack1996facs as "chromophore" that is the result of the F64L and S65T mutations from the #emph[Aequorea victoria] wild type GFP @prasher1992primary.

== Chromophore

Fluorescence in GFP begins with the absorption of a photon by the chromophore, typically in the blue region of the visible spectrum (around 488 nm).
This absorption process is intimately linked to the unique molecular structure of the chromophore.
The chromophore is formed autocatalytically from three amino acid residues---Ser65, Tyr66, and Gly67---through a series of reactions involving cyclization, dehydration, and oxidation.
The resulting structure consists of a hydroxybenzylidene imidazolinone moiety, which forms an extended $pi$-conjugated system.

#figure(
    image(FIG_DIR + "h-background/h005-cro/cro-a.svg", width: 2.3in),
    caption: [
        Neutral (i.e., A state) GFP chromophore.
    ],
    placement: auto
) <fig-cro-a>

=== Excitation

This π-conjugated system is crucial for the chromophore's light-absorbing properties.
The delocalized electrons in the conjugated bonds can be excited by photons of specific energies, corresponding to the energy gap between the ground state ($S_0$) and the first excited state ($S_1$) of the chromophore.
The exact absorption wavelength is fine-tuned by several factors that are discussed later: planarity, protonation state, and protein environment.

When a blue photon is absorbed, it promotes an electron from the highest occupied molecular orbital (HOMO) to the lowest unoccupied molecular orbital (LUMO) of the chromophore.
This electronic transition is predominantly $pi → pi^*$ in nature, reflecting the excitation within the π-conjugated system.

Following absorption, the excited chromophore undergoes rapid vibrational relaxation within the $S_1$ state, typically on a femtosecond to picosecond timescale.
This relaxation involves small structural adjustments in the chromophore and its immediate protein environment, preparing the system for the subsequent fluorescence emission.

=== Protonation states

The protonation state of the GFP chromophore is a critical determinant of its photophysical properties, playing a crucial role in the protein's spectral characteristics and fluorescence behavior.
The chromophore can exist in two primary forms: a neutral (protonated) state and an anionic (deprotonated) state, each exhibiting distinct spectroscopic signatures.

In its neutral form (shown in @fig-cro-a), the chromophore's phenolic oxygen is protonated, resulting in an absorption maximum typically around 395-400 nm.
This state is often referred to as the A state.
The neutral chromophore generally exhibits weaker fluorescence compared to its anionic counterpart, with emission maxima around 460 nm.
The reduced fluorescence efficiency of the neutral form is attributed to excited-state dynamics that favor non-radiative decay pathways.

The anionic form of the chromophore, where the phenolic oxygen is deprotonated, is primarily responsible for the characteristic green fluorescence of GFP.

#figure(
    image(FIG_DIR + "h-background/h005-cro/cro-b.svg", width: 2.3in),
    caption: [
        Anionic (i.e., B state) GFP chromophore.
    ],
    placement: auto
) <fig-cro-b>

This state, often called the B state, has an absorption maximum at approximately 475-490 nm and emits strongly at around 510 nm.
The anionic chromophore demonstrates a higher fluorescence quantum yield, making it the predominant contributor to GFP's bright fluorescence.

The relative population of these two states is influenced by several factors, including the local pH, specific interactions within the protein environment, and mutations in the protein sequence.
In wild-type GFP, the chromophore exists in an equilibrium between these two states, with the population distribution heavily dependent on pH.
Under physiological conditions, the anionic form is typically favored.

=== Planarity

The planarity of the eGFP chromophore plays a pivotal role in determining its fluorescent properties, particularly through enhanced conjugation and increased quantum yield.
These factors contribute significantly to the chromophore's spectroscopic characteristics and efficiency.

Enhanced conjugation in a planar chromophore structure is primarily due to the maximized overlap of p-orbitals in the π-conjugated system.
This optimal alignment of p-orbitals, perpendicular to the molecular plane, facilitates efficient delocalization of π-electrons across the entire conjugated system.
Such extensive electron delocalization has profound effects on the chromophore's electronic structure, most notably in reducing the energy gap between the highest occupied molecular orbital (HOMO) and the lowest unoccupied molecular orbital (LUMO).

Planarity significantly impacts the chromophore's quantum yield by restricting non-radiative decay pathways.
In non-planar configurations, rotation around single bonds, can serve as an efficient route for non-radiative decay.
However, when the chromophore is held in a planar conformation by the protein environment, these rotational movements are severely limited.
This restriction of molecular motion creates energy barriers in the excited state potential energy surface, effectively preventing the chromophore from accessing geometries that favor non-radiative decay.

Furthermore, the planar structure reduces coupling between electronic and vibrational states, which might otherwise lead to non-radiative relaxation through internal conversion.
The minimized structural changes between ground and excited states in a planar chromophore also contribute to reduced internal conversion rates.
Consequently, with fewer available non-radiative pathways, the excited chromophore is more likely to return to the ground state via fluorescence emission, directly increasing the fluorescence quantum yield.

=== Neutral-state fluorescence

Excitation and subsequent emission of the neutral state of the chromophore without any structural changes is one---infrequent---possibility.

Excitation of the neutral state occurs at approximately 395 nm, corresponding to the absorption of violet-blue light.
This excitation promotes the chromophore to its first excited singlet state ($S_1$) without immediate proton transfer.
The subsequent emission from this excited neutral state results in weak blue fluorescence with a peak around 460 nm.

#figure(
    image(FIG_DIR + "h-background/h006-cro-excitation/gfp-a-emission.png", width: 2.0in),
    caption: [
        Excitation and emission of neutral chromophore.
    ],
    placement: auto
) <fig-a-emission>

This protonated chromophore can be further stabilized by hydrogen bonding with a water molecule, Thr203, or Ser205.
The presence or absence of these interactions contribute to the distinct excitation and emission.
Although, this emission is relatively weak compared to the anionic state discussed next.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-a2-other.png", width: 3.0in),
    caption: [
        Example stabalizing configuration for the neutral chromophore.
    ],
    placement: auto
) <fig-a-residues>

=== Anionic-state fluorescence

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

The anionic state is stabilized by specific interactions within the protein barrel.
Notably, Thr203 plays a crucial role in stabilizing the anionic form through hydrogen bonding with the deprotonated phenolic oxygen.
Additionally, a protonated Glu222 could form a hydrogen bond with the anionic chromophore, further contributing to its stabilization.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-b-other.png", width: 3.0in),
    caption: [
        Example stabalizing configuration for the anionic chromophore.
    ],
    placement: auto
) <fig-b2-emission>

The protein environment around the chromophore is critical in maintaining this anionic configuration.
The β-barrel structure of eGFP provides a hydrophobic pocket that shields the chromophore from bulk solvent, contributing to the high quantum yield of fluorescence in this state.

=== Excited-state proton transfer

The excited-state proton transfer (ESPT) mechanism represents a fundamental process in GFP, contributing significantly to its unique spectroscopic properties.
This process involves the initial excitation of the neutral (protonated) chromophore, followed by rapid, successive proton transfer events in the excited state, ultimately resulting in emission from an anionic species.
An overview of the process is shown below.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-photocycle.png", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
) <fig-espt-emission>

1. Upon absorption of a photon at approximately 395 nm, the neutral chromophore is promoted to its first excited singlet state ($S_1$).

2. In this excited state, the chromophore exhibits markedly different acid-base properties compared to its ground state, becoming a much stronger acid.
    This enhanced acidity facilitates the transfer of a proton from the chromophore to a proximal acceptor within the protein matrix.
    The proton transfer pathway involves a sophisticated hydrogen-bonding network.
    A critical component of this network is a strategically positioned water molecule, which acts as the initial proton acceptor.
    This water molecule is part of a proton wire that includes Ser205 and terminates at Glu222. The transfer occurs on an ultrafast timescale, typically in the order of picoseconds.

3. Following the ESPT, the system exists transiently in an intermediate state (I\*), characterized by an anionic chromophore and a protonated Glu222.
    I\* subsequently relaxes and emits fluorescence at approximately 508 nm, closely resembling the emission profile of the intrinsically anionic chromophore.

4. The final step is a reverse ground-state proton transfer from Glu222 through Ser205, a water molecule, and terminated at the chromophore.


= roGFP2 contains redox-sensing cysteines

Redox-sensitive green fluorescent proteins (roGFPs) are engineered variants of GFP designed to report cellular redox states through changes in their fluorescence properties @hanson2004investigating.
Among these, roGFP2 has emerged as a particularly useful probe due to its ratiometric readout and midpoint potential, which are suitable for measuring redox conditions in reducing cellular compartments like the cytosol and mitochondria.

roGFP2 contains two cysteine residues (S147C and Q204C) on adjacent β-strands near the chromophore of a GFP variant already containing mutations C48S, S65T, and Q80R.
(A structural depiction of the relevant residues is in @fig-rogfp2-structure.)
The strategically placed cysteines can form a reversible disulfide bond in response to the surrounding redox environment changes.
Formation of this disulfide alters the protonation state of the chromophore, resulting in reciprocal changes in the excitation peaks at ~400 nm and ~490 nm.

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
We will refer to the neutral and anionic forms of the chromophore as "A state" and "B state", respectively.

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

Changes in A- or B-band absorption could indicate a variety of environmental changes in the chromophore.

- *Equilibrium ratio of A or B state populations.*
    Changes in neutral (A state) or anionic (B state) chromophore stability impacts the relative proportion of A- and B-band absorbance.
    For example, roGFP1 has a higher A band instead of B band @hanson2004investigating.
    A-state stability is directly correlated to A-band fluorescence; whereas the B band would be inversely correlated.
- *Excited-state proton transfer (ESPT) from A\* $arrow.r$ I\*.*
    Emissions at the typical 511 nm (green) fluorescence from the A state requires an ESPT from Cro66 to Glu222 through a coordinated water molecule and Ser205.
    Prohibiting ESPT would result in radiative emission at ~460 nm which often is not monitored experimentally.
- *Ground-state proton transfer (GSPT) from I $arrow.r$ A.*
    Reprotonating the chromophore through a GSPT is crucial for maintaining the A band and B-band lifetime.
    Disrupting the Glu222 $arrow.r$ Ser205, Ser205 $arrow.r$ H#sub[2]O, or H#sub[2]O $arrow.r$ Cro66 pathway would decrease the A-state population&mdash;likely with a corresponding B state increase.
- *Non-radiative emissions.*
    Enahnced flexibility of the chromophore through protein conformational shifts would lead to additional non-adiabatic crossings; thereby lowering the relative fluorescence in that state.

= Molecular simulations

== Protein preparation

Initial protein structures for reduced (#link("https://files.rcsb.org/download/1JC0.pdb")[1JC0]) and oxidized (#link("https://files.rcsb.org/download/1JC1.pdb")[1JC1]) states of roGFP2 were retrieved from the Protein Data Bank (PDB).
The structures were processed using in-house Python and bash scripts (available free of charge at #link("https://gitlab.com/oasci/studies/metalflare")[gitlab.com/oasci/studies/metalflare]).
The first chain, along with the crystallographic water molecules, were centered to the origin while minimize the box values using NumPy #addcite(), SciPy #addcite(), and MDAnalysis #addcite() packages.

Given the relevance and availability of the GFP mechanism's force field parameters, the chromophore was modeled in its anionic state.
All selenomethionine residues (MSE) were converted to methionine (MET), and Cys147 and Cys204 were transformed into the appropriate residues.
Glu222 was protonated to model the I state of the GFP photocycle; this still probes the anionic chromophore's stability while offering insight into GSPT dynamics.
The protonation states of all other residues were determined with PDB2PQR @jurrus2018improvements, using the default parameters to ensure standardization.
The pdb4amber tool was then used to validate the PDB file before proceeding.

== Simulation preparation

System preparation was performed using the tleap module of AmberTools v23.6 #addcite().
The protein structure was parameterized using the ff19SB force field #addcite().
For the solvent environment, we employed the OPC3 water model, #addcite(), which is known for its balanced representation of water properties in biomolecular simulations.
The 12-6 nonbonded model and parameters for all ions were taken from Sengupta et al @sengupta2021parameterization.
The system was neutralized by adding Na+ and Cl- ions as needed.
Additional ions were introduced to achieve a solvent ionic strength of 0.150 M to mimic physiological conditions.
The protein was solvated in a rectangular box, with a minimum distance of 10 Å between the protein and the box edges to minimize periodic boundary condition artifacts.
Parameters from Breyfogle et al. @breyfogle2023molecular were employed for the anionic chromophore.

== Minimization

The prepared system underwent a four-stage energy minimization protocol using Amber23 #addcite() to relieve any unfavorable interactions and optimize the structure.
All minimization stages used the steepest descent method for the first 1000 steps, followed by the conjugate gradient method for the remaining steps, with a maximum of 5000 steps per stage.
A non-bonded cutoff of 10.0 Å was applied throughout.
Periodic boundary conditions were employed, and coordinates were wrapped to the primary unit cell.
The minimization progress was monitored by writing energies every step and coordinates every 200 steps.

#emph[Stage 1:] Initial minimization was performed with restraints (force constant: 5.0 kcal/mol/Å²) on all non-hydrogen atoms of the entire system, allowing hydrogen atoms to relax and adjust their positions.
#emph[Stage 2:] The system was further minimized with restraints (force constant: 5.0 kcal/mol/Å²) on all non-hydrogen atoms of the solute (excluding water molecules and ions), allowing solvent and ions to equilibrate around the solute.
#emph[Stage 3:] Minimization continued with reduced restraints (force constant: 2.0 kcal/mol/Å²) applied only to the protein backbone, allowing side chains and other flexible parts to relax.
#emph[Stage 4:] Final minimization was performed with further reduced restraints (force constant: 1.0 kcal/mol/Å²).
The resulting minimized structure served as the starting point for subsequent relaxation and production simulations.

== Relaxation simulations

Following energy minimization, the system underwent a three-stage relaxation protocol using Amber23 to gradually equilibrate the structure and solvent.
Three independent runs were initiated with random initial velocities to ensure adequate sampling.
All subsequent simulations were continued using the respective run's restart files.

#emph[Stage 1:] An initial 20 ps NVT (constant Number of particles, Volume, and Temperature) simulation was performed with a 2 fs time step. The system was heated from 100 K to 300 K using Langevin dynamics with a collision frequency of 5 ps⁻¹. Restraints (force constant: 1.0 kcal/mol/Å²) were applied to the protein backbone.
SHAKE algorithm was used to constrain bonds involving hydrogen atoms.
The non-bonded cutoff was set to 10.0 Å.
#emph[Stage 2:] A 1 ns NPT simulation followed, maintaining the temperature at 300 K using Langevin dynamics (collision frequency: 5 ps⁻¹).
Pressure was regulated at 1.01325 bar using the Monte Carlo barostat with a relaxation time of 1 ps.
Restraints on the same atoms were reduced (force constant: 0.5 kcal/mol/Å²).
#emph[Stage 3:] The final relaxation stage consisted of a 1 ns NPT simulation with all positional restraints removed.

== Production simulations

All production runs were performed under the same setup as the last relaxation stage.
Each run was simulated for 500 ns with coordinates saved every ten ps.
The resulting trajectories from all three replicates were used for subsequent analyses, providing a cumulative 1.5 $mu$s of simulation data for the system under study.

== Hydrogen bond cutoff

A hydrogen bond of X---H $dots.c$ Y---Z, where X is the donor and Y is the acceptor atom, can be classified based on distances and angles.
One characteristic recommended by IUPAC is that the H $dots.c$ Y distance is less than the sum of H and Y van der Waals radii.
Hydrogen (1.10 Å) and oxygen (1.52 Å) @mantina2009consistent would have a cutoff of 2.62 Å.
Others @hubbard2010hydrogen recommend a cutoff of 2.50 Å based on structural analysis @mcdonald1994satisfying and quantum chemical calculations @liu2008geometrical.
Since the difference between a 2.5 and 2.62 Å cutoff is likely a substantially weak hydrogen bond, we will use a H $dots.c$ Y cutoff of 2.5 Å.

= Fluorescence mechanism of Cu(I) distinct from oxidation

Several experiments were performed to probe the fluorescence mechanism and binding affinity of Cu(I) to roGFP2.
Our focus here is to elucidate the distinct Cu(I) atomistic mechanism from the oxidized state.
@fig-cu-fluorescence shows the fluorescence of roGFP2 under various conditions.

#figure(
    image("s1b.png", width: 3.5in),
    caption: [
        Relative fluorescence (measured 528 nm emissions) of roGFP2 under reduced, oxidized, and Cu(I) conditions from 380 to 500 nm excitation scan.
        #todo("Check emission value?")
        Apo (i.e., reduced) roGFP2 exhibits typical bimodal absorption of A-band (excited at 400 nm) and B-band (excited at 488) peaks.
        Upon roGFP2 oxidation from 1 mM H#sub[2]O#sub[2], a shift in A- (increased) and B-band (decreased) absorption and subsequent 528 nm emission marks a corresponding change in neutral and anionic chromophore populations.
        Binding of Cu(I), however, exhibits a larger decrease in B-band without the A-band increase observed when oxidized.
    ],
    placement: auto
) <fig-cu-fluorescence>

First, we see a distinct difference from reduced $arrow.r$ oxidized and reduced $arrow.r$ Cu(I).
Oxidized roGFP2 results in a (1) enhanced A band and (2) slight decrease in the B band.
However, Cu(I) binding dramatically reduces the B band without change in the A band.
This indicates that Cu(I) binding affects the chromophore differently than oxidation.

= Cu(I) binding enhances roGFP2 backbone flexibility

First, we investigate the structural dynamics of Cys147 and Cys204 interactions by analyzing the C#sub[$alpha$]-C#sub[$alpha$] distances.
Experimental structures of both the reduced (PDB ID: #link("https://www.rcsb.org/structure/1JC0")[1JC0]) and oxidized (PDB ID: #link("https://www.rcsb.org/structure/1JC1")[1JC1]) states of for roGFP2 exhibited a mean C#sub[$alpha$]-C#sub[$alpha$] distance of 4.30 ± 0.12 and 4.07 ± 0.09, respectively @hanson2004investigating.
Our MD simulations agreed well with experimental observations as shown in Table 1.
Larger standard deviations, $sigma$, would be expected for molecular simulations.

#figure(
    table(
        columns: (auto, auto, auto),
        stroke: none,
        [*State*], [*Experimental (Å)*], [*MD simulations (Å)*],
        table.hline(),
        [Reduced], [4.30 ± 0.12], [4.34 ± 0.47],
        [Oxidized], [4.07 ± 0.09], [4.11 ± 0.29],
        [Cu(I)], [N/A], [4.78 ± 0.82],
    ),
    placement: auto
) <tab-alpha-c>

@fig-alpha-c shows the observed distribution of C#sub[$alpha$]-C#sub[$alpha$] distances.
Cu(I) binding to Cys147 and Cys204 in roGFP2 induces significant structural changes, particularly in the protein's conformation.
The observed increase in the C#sub[$alpha$]-C#sub[$alpha$] distance from approximately 4.3 Å to a broader distribution centered around 4.48 Å and 4.96 Å indicates a marked increase in conformational flexibility.

#figure(
    image(FIG_DIR + "b-cys/b004-cys147_ca-cys204_ca/b004-cys147_ca-cys204_ca-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
) <fig-alpha-c>

Cys147 is located near the C-terminus end of a $beta$-pleaded sheet between His148 and Thr203.
(Figure 7 illustrates the key residues in roGFP2 mechanism.)
We observe in Figure 8 that Cu(I) binding breaks this $beta$-sheet hydrogen bond.

In fact, Table 2 shows that the hydrogen bond probability ( &ndash;NH to O= within 2.5 Å) decreases from 0.865 (reduced) to 0.063 in Cu(I) simulations.

#figure(
    image(FIG_DIR + "f-backbone/f001-his148_h-thr203_o/f001-his148_h-thr203_o-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
)

#figure(
    caption: [Hydrogen bonding probability between residue backbones],
    table(
        columns: (auto, auto, auto),
        stroke: none,
        [*State*], [*His148 - Thr203*], [*Asn146 - Ser205*],
        table.hline(),
        [Reduced], [0.865], [0.047],
        [Oxidized], [0.997], [0.144],
        [Cu(I)], [0.063], [0.000],
    ),
    placement: auto
)

Parallel results are observed on the other side of Cys147&mdash;between Asn146 and Ser205.

#figure(
    image(FIG_DIR + "f-backbone/f002-ser205_h-asn146_o/f002-ser205_h-asn146_o-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
)

= Perturbations in anionic chromophore stability

#figure(
    image(FIG_DIR + "g-cro-interact/g006-cro66_oh-his148_hd1/g006-cro66_oh-his148_hd1-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
)

#figure(
    image(FIG_DIR + "g-cro-interact/g010-cro66_oh-thr203_hg1/g010-cro66_oh-thr203_hg1-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
)

= Disruption of ground-state proton transfer

#figure(
    image(FIG_DIR + "e-proton-wire/e001-ser205_og-glu222_he2/e001-ser205_og-glu222_he2-pdf.svg", width: 3.5in),
    caption: [
        #todo("Add caption")
    ],
    placement: auto
)
