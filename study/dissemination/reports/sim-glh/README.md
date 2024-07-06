# R002 - Molecular simulations with GLH222

We performed classical molecular dynamics (MD) simulations of roGFP2 with

-   the deprotonated chromophore in the ground state;
-   protonated GLU222 (abbreviated as GLH222).

These simulations allow us to probe the affect of roGFP2 oxidation and Cu(I) binding to CYS 147 and 204.
We model three states:

-   **Reduced**: CYS147 and CYS204 remain in their reduced form (i.e., no disulfide bridge) as our reference state;
-   **Oxidized**: CYS147 and CYS204 form a disulfide bond;
-   **Copper**: Cu(I) is coordinated to reduced CYS147 and CYS204.

Differences between **reduced** and **oxidized** states are indicative of the canonical roGFP2 atomistic mechanism.
If we observe similar differences between **reduced** and **copper** states, this would suggest a similar mechanism as the **oxidized** state; however, deviations between **oxidized** and **copper** would suggest a separate mechanism.

## Photocycle reprotonation

The last step of the [canonical photocycle](../fluorescence-mechanism/#photocycle) is reprotonation of of anionic chromophore through a GLU222, SER205, WAT, and finally CRO66 proton wire.
Intermolecular distances distances between each chemical species is a common feature to describe potential proton transfers in MD simulations.

### GLH222 to SER205

!!! quote "Figure X. Probability density of GLH222 to SER205 proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e008-ser205_og-glu222_he2/e008-ser205_og-glh222_he2-pdf.svg)
    </figure>

### SER205 to WAT

!!! quote "Figure X. Probability density of SER205 to WAT proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e010-ser205_og-h2o_h/e010-ser205_og-h2o_h-pdf.svg)
    </figure>

### WAT to CRO66

!!! quote "Figure X. Probability density of WAT to CRO66 proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e009-cro66_oh-h2o_h/e009-cro66_oh-h2o_h-pdf.svg)
    </figure>
