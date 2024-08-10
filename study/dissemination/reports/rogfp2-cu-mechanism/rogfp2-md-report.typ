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

Here, we propose atomistic mechanisms for a redox-sensitive green fluorescent protein (roGFP) variant including a novel Cu(I) sensing application.

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
    image(FIG_DIR + "h-background/h005-cro/cro-b.svg", width: 2.0in),
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
Such extensive electron delocalization has profound effects on the chromophore's electronic structure, most notably in reducing the energy gap between the Highest Occupied Molecular Orbital (HOMO) and the Lowest Unoccupied Molecular Orbital (LUMO).

The consequence of this reduced HOMO-LUMO gap is a red-shift in both absorption and emission spectra.
As the energy required for electronic transitions decreases, the chromophore absorbs and emits light at longer wavelengths.
This spectral shift is a direct result of the enhanced conjugation afforded by the planar structure.
Moreover, the planar configuration typically leads to higher oscillator strengths for electronic transitions, potentially resulting in stronger absorption and brighter fluorescence.

Planarity also significantly impacts the chromophore's quantum yield by restricting non-radiative decay pathways.
In non-planar configurations, rotation around single bonds, can serve as an efficient route for non-radiative decay.
However, when the chromophore is held in a planar conformation by the protein environment, these rotational movements are severely limited.
This restriction of molecular motion creates energy barriers in the excited state potential energy surface, effectively preventing the chromophore from accessing geometries that favor non-radiative decay.

Furthermore, the planar structure reduces coupling between electronic and vibrational states, which might otherwise lead to non-radiative relaxation through internal conversion.
The minimized structural changes between ground and excited states in a planar chromophore also contribute to reduced internal conversion rates.
Consequently, with fewer available non-radiative pathways, the excited chromophore is more likely to return to the ground state via fluorescence emission, directly increasing the fluorescence quantum yield.

=== Available mechanisms

With the above background, we will discuss possible fluorescence mechanisms.

==== Neutral state

Excitation and subsequent emission of the neutral state of the chromophore without any structural changes is one---infrequent---possibility.

Excitation of the neutral state occurs at approximately 395 nm, corresponding to the absorption of violet-blue light.
This excitation promotes the chromophore to its first excited singlet state ($S_1$) without immediate proton transfer.
The subsequent emission from this excited neutral state results in weak blue fluorescence with a peak around 460 nm.

#figure(
    image(FIG_DIR + "h-background/h006-cro-excitation/gfp-a-emission.svg", width: 2.0in),
    caption: [
        #todo("Add citation")
    ],
    placement: auto
) <fig-a-emission>

This protonated chromophore can be further stabilized by hydrogen bonding with a water molecule, Thr203, or Ser205.
The presence or absence of these interactions contribute to the distinct excitation and emission.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-a2.svg", width: 2.0in),
    caption: [
        #todo("Add citation")
    ],
    placement: auto
) <fig-a2-emission>

Although, this emission is relatively weak compared to the anionic state discussed next.

==== Anionic state

The anionic state of the chromophore represents a key configuration responsible for the protein's characteristic green fluorescence.
In this state, the chromophore exists in its deprotonated form, with the phenolic oxygen carrying a negative charge.

Excitation of the anionic chromophore occurs at approximately 475 nm, corresponding to the absorption of blue light.
This excitation promotes the chromophore to its first excited singlet state ($S_1$).
The subsequent relaxation and emission result in the bright green fluorescence typically associated with GFP, with a peak around 508 nm.

#figure(
    image(FIG_DIR + "h-background/h006-cro-excitation/gfp-b-emission.svg", width: 2.0in),
    caption: [
        #todo("Add citation")
    ],
    placement: auto
) <fig-b-emission>

The anionic state is stabilized by specific interactions within the protein barrel.
Notably, Thr203 plays a crucial role in stabilizing the anionic form through hydrogen bonding with the deprotonated phenolic oxygen.
Additionally, a protonated Glu222 could form a hydrogen bond with the anionic chromophore, further contributing to its stabilization.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-b.svg", width: 2.0in),
    caption: [
        #todo("Add citation")
    ],
    placement: auto
) <fig-b2-emission>

The protein environment around the chromophore is critical in maintaining this anionic configuration.
The β-barrel structure of eGFP provides a hydrophobic pocket that shields the chromophore from bulk solvent, contributing to the high quantum yield of fluorescence in this state.

==== Excited-state proton transfer

The excited-state proton transfer (ESPT) mechanism represents a fundamental process in GFP, contributing significantly to its unique spectroscopic properties.
This process involves the initial excitation of the neutral (protonated) chromophore, followed by rapid, successive proton transfer events in the excited state, ultimately resulting in emission from an anionic species.
An overview of the process is shown below.

#figure(
    image(FIG_DIR + "h-background/h004-cro-states/gfp-photocycle.svg", width: 2.0in),
    caption: [
        #todo("Add citation")
    ],
    placement: auto
) <fig-espt-emission>

1.  Upon absorption of a photon at approximately 395 nm, the neutral chromophore is promoted to its first excited singlet state ($S_1$).

2.  In this excited state, the chromophore exhibits markedly different acid-base properties compared to its ground state, becoming a much stronger acid.
    This enhanced acidity facilitates the transfer of a proton from the chromophore to a proximal acceptor within the protein matrix.

    The proton transfer pathway involves a sophisticated hydrogen-bonding network.
    A critical component of this network is a strategically positioned water molecule, which acts as the initial proton acceptor.
    This water molecule is part of a proton wire that includes Ser205 and terminates at Glu222. The transfer occurs on an ultrafast timescale, typically in the order of picoseconds.

3.  Following the ESPT, the system exists transiently in an intermediate state (I\*), characterized by an anionic chromophore and a protonated Glu222.
    I\* subsequently relaxes and emits fluorescence at approximately 508 nm, closely resembling the emission profile of the intrinsically anionic chromophore.

4.  The final step is a reverse ground-state proton transfer from Glu222 through Ser205, a water molecule, and terminated at the chromophore.

Little consensus has been made on the ESPT mechanism in GFP.










Most experimental studies characterizing GFP fluorescence will perform a scan of excitation wavelengths from 300 to 500 nm while monitoring around 511 nm.


GFP exhibits two distinct excitation peaks: (1) the A band around 400 nm and (2) the B band around 490 nm.
These peaks correspond to 511 nm fluorescence of the chromophore's excited neutral (protonated) and anionic (deprotonated) forms.
We will refer to the neutral and anionic forms of the chromophore as "A state" and "B state", respectively.

= roGFP2 contains redox-sensing cysteines

Redox-sensitive green fluorescent proteins (roGFPs) are engineered variants of GFP designed to report cellular redox states through changes in their fluorescence properties @hanson2004investigating.
Among these, roGFP2 has emerged as a particularly useful probe due to its ratiometric readout and midpoint potential, which are suitable for measuring redox conditions in reducing cellular compartments like the cytosol and mitochondria.

roGFP2 contains two cysteine residues (S147C and Q204C) on adjacent β-strands near the chromophore of a GFP variant already containing mutations C48S, S65T, and Q80R.
(A structural depiction of the relevant residues is in @fig-rogfp2-structure.)
The strategically placed cysteines can form a reversible disulfide bond in response to the surrounding redox environment changes.
Formation of this disulfide alters the protonation state of the chromophore, resulting in reciprocal changes in the excitation peaks at ~400 nm and ~490 nm.

#figure(
    image(FIG_DIR + "h-background/h009-rogfp2-sims/relevant-reduced.png", width: 5.0in),
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

= Methods

