# 001-rogfp-md

This experiment performs classical MD simulations of the [reduced form of roGFP](../../02-methods/01-protocols/gfp-definitions.md#reduced-form) for benchmarking structural changes of [mseGFP mutations](../../02-methods/01-protocols/gfp-definitions.md#mseGFP).

## Protein preparation

Protein Data Bank (PDB) files are not immediately usable for MD simulations.
Thus, we have to perform several steps to clean and prepare our PDB files.
Our protocol is located in [01-protein-prep.sh][protein-prep].

### Logging

First, we specify some environmental variables for metalflare.
This is mostly to enable logging.

```bash title="01-protein-prep.sh" linenums="7"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:7:13"
```

### Downloading

First, we download [roGFP](https://www.rcsb.org/structure/1JC0) from [RCSB](https://www.rcsb.org/).

```bash title="01-protein-prep.sh" linenums="20"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:20:21"
```

<div
    style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs'
    data-pdb='1jc0' data-backgroundalpha='0.0'
    data-select1='chain:A' data-style1='cartoon:color=spectrum'
    data-select2='resn:CRO' data-style2='cartoon:color=spectrum;stick'
    data-select3='resi:147' data-style3='cartoon:color=spectrum;stick'
    data-select4='resi:204' data-style4='cartoon:color=spectrum;stick'
    data-select5='chain:B' data-style5=''
    data-select6='chain:C' data-style6=''
    data-zoomto='chain:A'>
</div>

### Cleaning

```bash title="01-protein-prep.sh" linenums="23"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:23:28"
```

### pdb2pqr

```bash title="01-protein-prep.sh" linenums="30"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:30:30"
```

```bash title="01-protein-prep.sh" linenums="33"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:33:39"
```

```bash title="01-protein-prep.sh" linenums="40"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:40:40"
```

### pdb4amber

```bash title="01-protein-prep.sh" linenums="41"
--8<-- "study/03-data/001-rogfp-md/01-protein-prep.sh:42:42"
```

<!-- LINKS -->

[protein-prep]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/03-data/001-rogfp-md/01-protein-prep.sh

<!-- SCRIPTS -->

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
