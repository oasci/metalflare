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

= Methods

== Molecular simulations

Four simulation sets---reduced, oxidized, Cu(I), and Na#super[+]---were performed to gain atomistic insight into roGFP2 sensing mechanisms.
The reduced simulations represent our designed control system with protonated Cys147 and Cys204.
We will refer to this crucial cysteine pair as the 'sensing cysteines'.
Oxidized simulations contain a disulfide bond linking the sensing cysteines.
Cu(I) simulations let the ion coordinate with the sensing cysteines and use their respective metal-coordinated parameters.
Na#super[+] simulations use the same Cys metal-coordinated parameters without adding additional ions, which offers a control for Cys parameters and ion identity.

== Protein preparation

Initial protein structures for reduced (#link("https://files.rcsb.org/download/1JC0.pdb")[1JC0]) and oxidized (#link("https://files.rcsb.org/download/1JC1.pdb")[1JC1]) states of roGFP2 #addcite() were retrieved from the Protein Data Bank (PDB) #addcite().
The structures were processed using in-house Python and bash scripts (available free of charge at #link("github.com/oasci/metalflare")[github.com/oasci/metalflare]).
The first chain, along with the crystallographic water molecules, were centered to the origin and rotated to minimize the box size using NumPy #addcite(), SciPy #addcite(), and MDAnalysis #addcite() packages.

Given the relevance and availability of the GFP mechanism's force field parameters, the chromophore was modeled in its anionic state.
All selenomethionine residues (MSE) were converted to methionine (MET), and Cys147 and Cys204 were transformed into the appropriate states.
Glu222 was protonated to model the I state of the GFP photocycle; this still probes the anionic chromophore's stability while offering insight into ground-state proton transfer (GSPT) dynamics.
The protonation states of all other residues were determined with PDB2PQR @jurrus2018improvements, using the default parameters.
pdb4amber tool #addcite() was then used to validate the PDB file before proceeding.

== Simulation preparation

System preparation was performed using the tleap module of AmberTools v23.6 #addcite().
The protein structure was parameterized using the ff19SB force field #addcite().
For the solvent environment, we employed the OPC3 water model, #addcite(), known for its balanced representation of water properties in biomolecular simulations.
The 12-6 nonbonded model and parameters for all ions were taken from Sengupta et al. @sengupta2021parameterization
The system was neutralized by adding Na#super[+] and Cl#super[-] ions as needed.
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

== Relaxation

Following energy minimization, the system underwent a three-stage relaxation protocol using Amber23 to equilibrate the structure and solvent gradually.
Three independent runs were initiated with random initial velocities drawn from the Maxwell-Boltzmann distribution at 100 K. #todo("Double check this")
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
The resulting trajectories from all three replicates were used for subsequent analyses, providing a cumulative 1.5 $mu$s of simulation data for each system state.


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
) <fig-1>
