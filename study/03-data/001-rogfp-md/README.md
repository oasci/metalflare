# 001-rogfp-md

MD simulations of the reduction-oxidation-sensitive GFP.

## Protein preparation

Protein Data Bank (PDB) files are not immediately usable for MD simulations.
Thus, we have to perform several steps to clean and prepare our PDB files.
All of these steps are located in the [protein-prep.sh script][protein-prep].

### Logging

First, we specify some environmental variables for metalflare.
This is mostly to enable logging.

```bash title="protein-prep.sh" linenums="7"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:7:12"
```

### Downloading

First, we download [roGFP](https://www.rcsb.org/structure/1JC0) from [RCSB](https://www.rcsb.org/).

```bash title="protein-prep.sh" linenums="19"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:19:20"
```

<div style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs' data-pdb='1JC0' data-backgroundalpha='0.0' data-select1='chain:A' data-style1='cartoon:color=spectrum' data-select2='resn:CRO' data-style2='stick' data-select3='chain:B' data-style3='hidden' data-select4='chain:C' data-style4='hidden'  data-zoomto='chain:A'></div>

### Cleaning

```bash title="protein-prep.sh" linenums="22"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:22:27"
```

### pdb2pqr

```bash title="protein-prep.sh" linenums="29"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:29:29"
```

```bash title="protein-prep.sh" linenums="32"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:32:38"
```

```bash title="protein-prep.sh" linenums="39"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:39:39"
```

### pdb4amber

```bash title="protein-prep.sh" linenums="41"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:41:41"
```

<div style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs' data-href='https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/03-data/001-rogfp-md/structures/protein/1JC0-final.pdb' data-backgroundalpha='0.0' data-select1='chain:A' data-style1='cartoon:color=spectrum' data-select2='resn:CRO' data-style2='stick'></div>

[protein-prep]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/03-data/001-rogfp-md/protein-prep.sh

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
