# R002 - Molecular simulations

We performed classical molecular dynamics (MD) simulations of roGFP2 [^hanson2004investigating] with

-   the anionic chromophore in the ground state;
-   neutral and anionic Glu222.

roGFP2 is derived from enhanced GFP (eGFP; S65T and F64L mutations) with S147C and Q204C mutations.
These simulations allow us to probe the affect of roGFP2 oxidation and Cu(I) binding of Cys147 and 204.
We performed the following three sets of simulations.

=== "Reduced"

    Cys147 and Cys204 remain in their reduced (i.e., protonated) form.

    <figure markdown>
    ![](../../../figures/h-background/h008-cys-aligned/cys-sensor-reduced.png){ width=500 }
    </figure>

=== "Oxidized"

    Disulfide bridge between Cys147 and Cys204.

    <figure markdown>
    ![](../../../figures/h-background/h008-cys-aligned/cys-sensor-oxidized.png){ width=500 }
    </figure>

=== "Copper"

    Reduced Cys147 and Cys204 with a coordinated Cu(I).

    <figure markdown>
    ![](../../../figures/h-background/h008-cys-aligned/cys-sensor-cu.png){ width=500 }
    </figure>

Differences between reduced and oxidized states are indicative of the canonical roGFP2 atomistic mechanism.
If we observe similar differences between reduced and copper states, this would suggest a similar mechanism as the oxidized state; however, deviations between oxidized and copper would suggest a separate mechanism.

## Cysteine conformations

We first investigate the structural dynamics of Cys147 and Cys204 interactions by analyzing the C$_\alpha$-C$_\alpha$ distances.

!!! quote ""
    <figure markdown>
    ![](../../../figures/f-cys-beta/f008h-cys147_ca-cys204_ca/f008h-cys147_ca-cys204_ca-pdf.svg){ width=600 }
    </figure>

    For more figure information, go [here](../../../figures/f-cys-beta/f008h-cys147_ca-cys204_ca/).

As one would expect, formation of the disulfide bridge in the oxidized state induces a highly strained conformation of around 4.07 Å C$_\alpha$-C$_\alpha$ distance, agreeing with the average experimental distance of 4.04 ± 0.09 Å [^hanson2004investigating].
Additional conformational flexibility is observed in the reduced simulations with a broader distribution centered around 5.07 Å.
A reduced&mdash;albeit partially oxidized&mdash;crystal structure of roGFP2 had a distance of 4.30 ± 0.12 Å [^hanson2004investigating], well within the observed distribution.
Our Cu(I) simulation C$_\alpha$-C$_\alpha$ distances resemble those of the reduced cysteines; however, Cu(I) binding induces closer contacts on the order of 0.5 Å.
Certainly, binding of Cu(I) to Cys147 and Cys204 causes a distinct conformational change relative to canonical roGFP2.

## Photocycle reprotonation

The last step of the [canonical photocycle](../fluorescence-mechanism/#photocycle) is reprotonation of of anionic chromophore through a GLU222, SER205, WAT, and finally CRO66 proton wire.
Intermolecular distances distances between each chemical species is a common feature to describe potential proton transfers in MD simulations.

### GLH222 to SER205

!!! quote "Figure X. Illustration of GLH222 to SER205 distance."
    <figure markdown>
    ![](../../../figures/h-background/h007-distances/gfp-glh222-ser205.svg){ width=500 }
    </figure>

!!! quote "Figure X. Probability density of GLH222 to SER205 proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e008-ser205_og-glu222_he2/e008-ser205_og-glh222_he2-pdf.svg){ width=500 }
    </figure>

### SER205 to WAT

!!! quote "Figure X. Illustration of SER205 to WAT distance."
    <figure markdown>
    ![](../../../figures/h-background/h007-distances/gfp-ser205-wat.svg){ width=500 }
    </figure>

!!! quote "Figure X. Probability density of SER205 to WAT proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e011-ser205_hg-h2o_o/e010-ser205_hg-h2o_o-pdf.svg){ width=500 }
    </figure>

### WAT to CRO66

!!! quote "Figure X. Illustration of WAT to CRO distance."
    <figure markdown>
    ![](../../../figures/h-background/h007-distances/gfp-wat-cro.svg){ width=500 }
    </figure>

!!! quote "Figure X. Probability density of WAT to CRO66 proton transfer."
    <figure markdown>
    ![](../../../figures/e-proton-wire/e009-cro66_oh-h2o_h/e009-cro66_oh-h2o_h-pdf.svg){ width=500 }
    </figure>

## B-state stabilization

TODO:

!!! quote "Figure X. Probability density of THR203 to CRO66 hydrogen bonding in GLH222 simulations."
    <figure markdown>
    ![](../../../figures/h-background/h007-distances/gfp-b-thr203-cro66.svg){ width=500 }
    </figure>

!!! quote "Figure X. Probability density of THR203 to CRO66 hydrogen bonding in GLH222 simulations."
    <figure markdown>
    ![](../../../figures/b-cro-between/b008h-cro66_oh-thr203_hg1/b008h-cro66_oh-thr203_hg1-pdf.svg){ width=500 }
    </figure>

### Deprotonated GLU

!!! quote "Figure X. Probability density of THR203 to CRO66 hydrogen bonding in GLU222 simulations."
    <figure markdown>
    ![](../../../figures/b-cro-between/b008-cro66_oh-thr203_hg1/b008-cro66_oh-thr203_hg1-pdf.svg){ width=500 }
    </figure>

## Non-adiabatic crossings

!!! quote "Figure X."
    <figure markdown>
    ![](../../../figures/h-background/h005-cro/cro-b-atom-types.svg){ width=500 }
    </figure>

!!! quote "Figure X."
    <figure markdown>
    ![](../../../figures/a-cro/a003-cro66_cd2_cg2_cb2_ca2/a003-cro66_cd2_cg2_cb2_ca2-pdf.svg){ width=500 }
    </figure>

!!! quote "Figure X."
    <figure markdown>
    ![](../../../figures/a-cro/a004-cro66_cg2_cb2_ca2_c2/a004-cro66_cg2_cb2_ca2_c2-pdf.svg){ width=500 }
    </figure>

<!-- References -->

[^hanson2004investigating]: Hanson, G. T., Aggeler, R., Oglesbee, D., Cannon, M., Capaldi, R. A., Tsien, R. Y., & Remington, S. J. (2004). Investigating mitochondrial redox potential with redox-sensitive green fluorescent protein indicators. Journal of Biological Chemistry, 279(13), 13044-13053. DOI: [10.1074/jbc.M312846200](https://doi.org/10.1074/jbc.M312846200)
