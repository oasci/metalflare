# GFP definitions

## eGFP

[Enhanced GFP][2y0g] (eGFP), first introduced by [Heim et al.][egfp paper], has `S65T` and `F64L` mutations from the wild type protein.

<div
    style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs'
    data-pdb='2Y0G' data-backgroundalpha='0.0'
    data-select1='chain:A' data-style1='cartoon:color=spectrum'
    data-select2='resn:CRO' data-style2='cartoon:color=spectrum;stick'
    data-select3='resi:147' data-style3='cartoon:color=spectrum;stick'
    data-select4='resi:204' data-style4='cartoon:color=spectrum;stick'
    data-zoomto='chain:A'>
</div>

## roGFP

The [redox-sensitive GFP][rogfp paper] (roGFP) is derived from [eGFP](#egfp) with two additional mutations: `S147C` and `Q204C`.
Introduced as roGFP2 [Hanson et al.][rogfp paper], this forms a reversible formation of a [reduced][1jc0] and [oxidized][1jc1] disulfide bridge between `147` and `204`.

### Reduced form

[1JC0][1jc0] shows the reduced (i.e., broken) form of `147`-`204` disulfide bond.

<div
    style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs'
    data-pdb='1JC0' data-backgroundalpha='0.0'
    data-select1='chain:A' data-style1='cartoon:color=spectrum'
    data-select2='resn:CRO' data-style2='cartoon:color=spectrum;stick'
    data-select3='resi:147' data-style3='cartoon:color=spectrum;stick'
    data-select4='resi:204' data-style4='cartoon:color=spectrum;stick'
    data-select5='chain:B' data-style5=''
    data-select6='chain:C' data-style6=''
    data-zoomto='chain:A;resi:204'>
</div>

### Oxidized form

[1JC1][1jc1] shows the oxidized (i.e., formed) form of `147`-`204` disulfide bond.

<div
    style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs'
    data-pdb='1JC1' data-backgroundalpha='0.0'
    data-select1='chain:A' data-style1='cartoon:color=spectrum'
    data-select2='resn:CRO' data-style2='cartoon:color=spectrum;stick'
    data-select3='resi:147' data-style3='cartoon:color=spectrum;stick'
    data-select4='resi:204' data-style4='cartoon:color=spectrum;stick'
    data-select5='chain:B' data-style5=''
    data-select6='chain:C' data-style6=''
    data-zoomto='chain:A;resi:204'>
</div>

## msGFP

Metal-sensing GFP (msGFP) is similar to roGFP, but mutates `S147C` and `S202C` from [eGFP][2y0g].
Since the protein can no longer make the `147`-`204` disulfide bond, all msGFP simulations will start from the [reduced form](#reduced-form) with `C202S` and `S202C` mutations.

<!-- LINKS -->

[egfp paper]: https://doi.org/10.1038/373663b0
[1jc0]: https://www.rcsb.org/structure/1jc0
[1jc1]: https://www.rcsb.org/structure/1jc1
[rogfp paper]: https://doi.org/10.1074/jbc.M312846200
[2y0g]: https://www.rcsb.org/structure/2y0g

<!-- SCRIPTS -->

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
