# 02 - Protein preparation from RCSB

Protein Data Bank (PDB) files are not immediately usable for MD simulations.
Thus, we have to perform several steps to clean and prepare our PDB files.
Our protocol is located in [01-protein-prep.sh][protein-prep.sh].

## Logging

First, we specify some environmental variables for metalflare.
This is mostly to enable logging.

```bash
export METALFLARE_LOG=True
export METALFLARE_STDOUT=True
export METALFLARE_LOG_LEVEL=10
```

## Download PDB

First, we download the protein from [RCSB](https://www.rcsb.org/).
For an example, we will be using [1JC0](https://www.rcsb.org/structure/1JC0).

```bash
wget https://files.rcsb.org/download/$PDB_ID.pdb -O $SAVE_DIR/0-$PDB_ID.pdb
```

We receive the following structure.

<div id="original-pdb" class="mol-container"></div>
<script>
var viewer1 = $3Dmol.createViewer(
    document.querySelector('#original-pdb'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/1JC0.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer1.addModel( data, 'pdb' );
        viewer1.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer1.setStyle({resn: 'CRO'}, {stick: {}});
        viewer1.setStyle({resn: 'HOH'}, {sphere: {scale: '0.3', opacity: '0.95'}});
        viewer1.zoomTo();
        viewer1.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

## Select protein

We will only be working with one protein and need to select chain `A`.
All water molecules (residue name `HOH`) are also removed.

We use the [`metalflare-select-atoms`][pdb.select.cli_select_atoms] script which just drives [`pdb.select.run_select_atoms()`][pdb.select.run_select_atoms].


```bash
metalflare-select-atoms $SAVE_DIR/0-$PDB_ID.pdb $SAVE_DIR/1-$PDB_ID-chain-A.pdb --select_str chainID A and not resname HOH
```
<div id="select-chain-a" class="mol-container"></div>
<script>
var viewer2 = $3Dmol.createViewer(
    document.querySelector('#select-chain-a'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/1JC0.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer2.addModel( data, 'pdb' );
        viewer2.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer2.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer2.setStyle({chain: 'B'}, {});
        viewer2.setStyle({chain: 'C'}, {});
        viewer2.zoomTo({chain: 'A'});
        viewer2.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

## Keep only `ATOM` and `HETATM` lines



## pdb2pqr

```bash title="01-protein-prep.sh" linenums="30"
--8<-- "study/03-data/001-rogfp-md/scripts/01-protein-prep.sh:30:30"
```

```bash title="01-protein-prep.sh" linenums="33"
--8<-- "study/03-data/001-rogfp-md/scripts/01-protein-prep.sh:33:39"
```

```bash title="01-protein-prep.sh" linenums="40"
--8<-- "study/03-data/001-rogfp-md/scripts/01-protein-prep.sh:40:40"
```

## pdb4amber

```bash title="01-protein-prep.sh" linenums="41"
--8<-- "study/03-data/001-rogfp-md/scripts/01-protein-prep.sh:42:42"
```

<!-- LINKS -->

[protein-prep.sh]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/03-data/001-rogfp-md/scripts/01-protein-prep.sh
