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
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:7:13"
```

### Downloading

First, we download [roGFP](https://www.rcsb.org/structure/1JC0) from [RCSB](https://www.rcsb.org/).

```bash title="protein-prep.sh" linenums="20"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:20:21"
```

<div style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs' data-pdb='1JC0' data-backgroundalpha='0.0' data-select1='chain:A' data-style1='cartoon:color=spectrum' data-select2='resn:CRO' data-style2='stick' data-select3='chain:B' data-style3='hidden' data-select4='chain:C' data-style4='hidden'  data-zoomto='chain:A'></div>

### Cleaning

```bash title="protein-prep.sh" linenums="23"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:23:28"
```

### pdb2pqr

```bash title="protein-prep.sh" linenums="30"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:30:30"
```

```bash title="protein-prep.sh" linenums="33"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:33:39"
```

```bash title="protein-prep.sh" linenums="40"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:40:40"
```

### pdb4amber

```bash title="protein-prep.sh" linenums="41"
--8<-- "study/03-data/001-rogfp-md/protein-prep.sh:42:42"
```

<div style="height: 400px; width: 100%; position: relative;" class='viewer_3Dmoljs' data-href='https://gitlab.com/oasci/studies/metalflare/-/raw/main/study/03-data/001-rogfp-md/structures/protein/1JC0-final.pdb' data-backgroundalpha='0.0' data-select1='chain:A' data-style1='cartoon:color=spectrum' data-select2='resn:CRO' data-style2='stick'></div>

[protein-prep]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/03-data/001-rogfp-md/protein-prep.sh

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
