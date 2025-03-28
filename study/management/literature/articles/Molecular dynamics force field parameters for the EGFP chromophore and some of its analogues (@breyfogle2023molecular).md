---
tags:
  - lit
date: 2023-06-26
---

# Molecular dynamics force field parameters for the EGFP chromophore and some of its analogues

**Authors:** Kimberly L. Breyfogle, Dalton L. Blood, Andreana M. Rosnik, Brent P. Krueger

**DOI:** [10.1021/acs.jpcb.3c01486](https://doi.org/10.1021/acs.jpcb.3c01486)

<!-- more -->

Takeaways:

-   Parameterizes chromophores for GFP.

## Abstract

> Fluorescent proteins (FPs) have had an enormous impact on molecular and cellular biology and are employed in a wide range of studies of molecular structure and dynamics. Yet, only a modest number of papers have published molecular dynamics (MD) parameters describing FPs. And despite the development of a wide range of FPs, there has been no careful development of MD parameters across a series of FPs. In this work, we present MD parameters describing six fluorescent protein chromophores (EGFP, EBFP, EYFP, ECFP, mCherry, and DsRed) for use with the Cornell _et al_. ( _J. Am. Chem. Soc._ 1995, _117_, 5179−5197) family of AMBER force fields, including ff14SB and ff19SB. We explore a wide range of solvent dielectric constants for determining the chromophore equilibrium geometry and evaluate the impact of the modeled solvent on the final atomic charges. We also present our methodological approach in which we considered all six chromophores together with a focus on modularity, transferability, and balance with existing force fields. The parameters given here make it easy to employ MD simulations to study any of the six systems, whereas the methodology makes it easy for anyone to extend this work to develop consistent parameters for additional fluorescent proteins. The results of our own MD simulations are presented, showing that the classical MD parameters yield chromophore structural distributions that compare well with QM/MM simulations.

## Introduction

> In this work, we report MD parameters for six different fluorescent protein chromophores (EGFP, EBFP, EYFP, ECFP, mCherry, and DsRed with two forms of EBFP; see Figure 2 for use with ff14SB (44) and ff19SB, (45) the two newest versions of the family of Cornell-based (21) AMBER force fields. (46)

## Methods

The eGFP chromophore (PDB code: `CRO`) used the PDB reference of `2Y0G` with a chromophore charge of `-1`.

> The leading and trailing residues were converted into acetyl (ACE) and _N_-methyl amide (NME) capping groups and the structures were then optimized using a four-step process.

> A wide range of solvent dielectrics (from 1.1 to 30) have been used to model protein environments in various studies over the past three decades, (67−80) although most of these works have utilized ε = 10–20.
> Because of this ambiguity, we tested the full optimization procedure described above using six different solvents, including
>
> -   vacuum (ε = 1),
> -   pentylamine (ε = 4.20),s
> -   1-octanol (ε = 9.86),
> -   1-pentanol (ε = 15.13),
> -   2-propanol (ε = 19.26), and
> -   water (ε = 78.36),
>  
> all implemented with the SMD solvent model. (55)
> These solvents were chosen based on a number of factors:
> 1.  in the crystal structure, there are a variety of sidechains within 3 Å of, and pointing toward, the chromophore including hydrophobic valine, aromatic/polar histidine, and multiple polar and/or charged groups including two threonines, arginine, glutamine, glutamic acid, and three waters, such that short-chain alcohols or amines are a reasonable compromise for mimicking this environment;
> 2.  surface tensions in the range of 30–37 cal/mol·Å2, consistent with alcohols and amines of modest chain length; and
> 3.  dielectric constants that span the range of 4–20 plus the two extreme values of vacuum and water.

> To maintain consistency and balance with the Cornell-based force fields, we sought to define the fewest possible number of new parameters. For those new parameters that would be needed, we chose to base them on either ff14/19SB or the generalized AMBER force field (gaff), (82) which is designed to be compatible with the Cornell family of force fields.

> Once atom types were determined, most bond, angle, dihedral, and van der Waals parameters were already defined by ff14SB or gaff and were used without adjustment. Any remaining (undefined) angle and dihedral parameters were assigned by analogy to ff14SB or gaff as described below in Results and Discussion.

> Atomic charges were determined using two-stage Restrained Electrostatic Potential (RESP) fits of electrostatic potentials (ESPs) that were computed on B3LYP/6-31G(d)/SMD(1-pentanol) optimized structures using HF/6-31G(d), as prescribed by Kollman and co-workers.

> The RESP module within the AMBER suite (47) was used to perform the fitting. Within the RESP method, we imposed three requirements consistent with the original Cornell charge model:
> -   the total charge of all capping group atoms must be zero to allow them to be removed while maintaining an integer charge on the portion of the chromophore that will serve as the residue;
> -   atomic charges of rotationally equivalent atoms within the molecule must be equal, with the exception of those found on capping groups; and
> -   the atomic charges of the terminal/amide N, H, C, and O atoms in each chromophore must be consistent with the corresponding atoms in the standard ff19SB amino acids.
