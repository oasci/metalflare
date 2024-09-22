# R003 - Cu(I) binding mechanism of roGFP2

<iframe width="100%" height="800" src="./rogfp2-cu-comp.pdf">

-   **Inadequate Modeling of Cu(I):** The Cu(I) ion presents unique challenges due to its d<sup>10</sup> electronic configuration and propensity for variable coordination geometries and soft interactions, which are not well-captured by classical force fields.
    The manuscript lacks a detailed justification for the chosen parameters from Sengupta et al. and does not discuss their suitability for modeling Cu(I) coordination with cysteine residues.
-   **Absence of Polarization Effects:** The simulations do not account for electronic polarization, which is critical in metal-ligand interactions.
    The use of non-polarizable force fields may lead to inaccurate representations of the Cu(I)-protein interactions, potentially compromising the validity of the results.
-   **Parameter Validation:** The parameters for the anionic chromophore are taken from Breyfogle et al., but there is no discussion on their validation or applicability to the current system,.
-   **Excited-State Considerations:** The study mentions limitations due to the lack of parameters for the excited and neutral chromophore states.
    This is a critical issue since fluorescence properties are inherently linked to the electronic excited states.
    Ignoring these states limits the ability to draw meaningful conclusions about fluorescence mechanisms.
-   **Simulation Length and Convergence:** While each production run is 500 ns with three replicates, it's unclear whether this duration is sufficient for the systems to reach equilibrium, especially for capturing slow conformational changes associated with metal binding and β-strand fraying.
-   **Simplistic Criteria:** The hydrogen bond definition relies solely on a distance cutoff of 2.5 Å, neglecting the angular dependence that is crucial for accurate hydrogen bonding characterization.
    This omission may lead to an overestimation of hydrogen bond populations and misinterpretation of the interactions.
-   **Limitations of Classical MD:** Classical MD simulations with fixed protonation states are incapable of modeling proton transfer events, which involve bond breaking and formation. The study's conclusions about GSPT disruption are speculative without quantum mechanical (QM) or reactive MD simulations.
-   **Correlation with Fluorescence Quenching:** The direct connection between the observed structural changes and fluorescence quenching is not convincingly established, especially given the absence of electronic excited-state modeling.
