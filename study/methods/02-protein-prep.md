# 02 - Protein preparation from RCSB

Protein Data Bank (PDB) files are not immediately usable for MD simulations.
Thus, we have to perform several steps to clean and prepare our PDB files.

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
var uri = 'https://files.rcsb.org/view/1JC0.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#original-pdb'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({resn: 'CRO'}, {stick: {}});
        viewer.setStyle({resn: 'HOH'}, {sphere: {scale: '0.3', opacity: '0.95'}});
        viewer.zoomTo();
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
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
var uri = 'https://files.rcsb.org/view/1JC0.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#select-chain-a'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.zoomTo({chain: 'A'});
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

## Minify PDB lines

Keep only `ATOM` and `HETATM` lines.

```bash
metalflare-filter-pdb $SAVE_DIR/1-$PDB_ID-chain-A.pdb  --output $SAVE_DIR/2-$PDB_ID-filtered.pdb
```

## Centering

```text
ATOM    471  CG  LEU A  64     172.089   9.780  36.151  1.00 35.35      A    C
ATOM    472  CD1 LEU A  64     172.939   9.531  34.924  1.00 28.06      A    C
ATOM    473  CD2 LEU A  64     170.620   9.387  35.892  1.00 33.01      A    C
HETATM  474  N1  CRO A  65     173.570   8.493  40.293  1.00 24.05      A    N
HETATM  475  CA1 CRO A  65     174.025   7.483  41.259  1.00 25.10      A    C
HETATM  476  CB1 CRO A  65     175.255   6.742  40.591  1.00 31.02      A    C
```

```bash
metalflare-center $SAVE_DIR/2-$PDB_ID-filtered.pdb --output $SAVE_DIR/3-$PDB_ID-centered.pdb
```

```text
ATOM    471  CG  LEU A  64      -1.838  -0.100  -6.847  1.00 35.35      A    C
ATOM    472  CD1 LEU A  64      -0.988  -0.349  -8.074  1.00 28.06      A    C
ATOM    473  CD2 LEU A  64      -3.307  -0.493  -7.106  1.00 33.01      A    C
HETATM  474  N1  CRO A  65      -0.357  -1.387  -2.705  1.00 24.05      A    N
HETATM  475  CA1 CRO A  65       0.098  -2.397  -1.739  1.00 25.10      A    C
HETATM  476  CB1 CRO A  65       1.328  -3.138  -2.407  1.00 31.02      A    C
```

## Rotate structure

Will attempt to rotate the structure in order to minimize the number of water molecules we will add later.

```bash
metalflare-minimize-box $SAVE_DIR/3-$PDB_ID-centered.pdb --output $SAVE_DIR/4-$PDB_ID-rotated.pdb
```

For our 1JC0 example, the box volume decreased from 88 063 to 71 986 Å<sup>3</sup> after this script, which decreases the number of water molecules by at least 500.

## Unify residue IDs

```text
ATOM    471  CG  LEU A  64     172.089   9.780  36.151  1.00 35.35      A    C
ATOM    472  CD1 LEU A  64     172.939   9.531  34.924  1.00 28.06      A    C
ATOM    473  CD2 LEU A  64     170.620   9.387  35.892  1.00 33.01      A    C
HETATM  474  N1  CRO A  66     173.570   8.493  40.293  1.00 24.05      A    N
HETATM  475  CA1 CRO A  66     174.025   7.483  41.259  1.00 25.10      A    C
HETATM  476  CB1 CRO A  66     175.255   6.742  40.591  1.00 31.02      A    C
```

```bash
metalflare-unify-resids $SAVE_DIR/4-$PDB_ID-rotated.pdb --output $SAVE_DIR/5-$PDB_ID-residues.pdb
```

```text
ATOM    471  CG  LEU A  64     172.089   9.780  36.151  1.00 35.35      A    C
ATOM    472  CD1 LEU A  64     172.939   9.531  34.924  1.00 28.06      A    C
ATOM    473  CD2 LEU A  64     170.620   9.387  35.892  1.00 33.01      A    C
HETATM  474  N1  CRO A  65     173.570   8.493  40.293  1.00 24.05      A    N
HETATM  475  CA1 CRO A  65     174.025   7.483  41.259  1.00 25.10      A    C
HETATM  476  CB1 CRO A  65     175.255   6.742  40.591  1.00 31.02      A    C
```

## Residue states

### Methionine

Methionine (`MET`) residues are often artificially changed to selenomethionine (`MSE`) to ensure proper crystallization by multi-wavelength anomalous dispersion.
We almost always want to model with the wild-type `MET`, so we replace any `MSE` with `MET` residues and `Se` atoms to `S`.

```bash
metalflare-rename-resname $SAVE_DIR/5-$PDB_ID-residues.pdb MSE MET --output $SAVE_DIR/5-$PDB_ID-residues.pdb
```

### Cysteine

<!-- Amber uses `CYX` instead of `CYS` to indicate that cysteine residues are involved in disulfide bonds.
Often this has to be manually inspected and changed.
If you want to convert all `CYX` to `CYS` residues to ensure no disulfide bonds are present, you can use the following script. -->

```bash
metalflare-rename-resname $SAVE_DIR/5-$PDB_ID-residues.pdb CYS CYM --include 145 202 --output $SAVE_DIR/5-$PDB_ID-residues.pdb
```

## Protonation and steric clashes

[PDB2PQR][pdb2pqr] predicts protonation states of histidine (`HIS`), aspartic acid (`ASP`), glutamic acid (`GLU`), lysine (`LYS`).

```bash
pdb2pqr --log-level INFO --ff=AMBER --keep-chain --ffout=AMBER $SAVE_DIR/5-$PDB_ID-residues.pdb $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
```

Sometimes [PDB2PQR][pdb2pqr] cannot process some atoms, so we need to add them back.

```bash
metalflare-merge-pdbs $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb $SAVE_DIR/5-$PDB_ID-residues.pdb --output $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb
```

!!! warning

    [PDB2PQR][pdb2pqr] cannot process non-standard residues (e.g., the GFP chromophore) and thus cannot add hydrogens to them.
    These are often added later using a program like tleap.

## Unify water residues

```bash
metalflare-rename-resname $SAVE_DIR/6-$PDB_ID-pdb2pqr.pdb HOH WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb
metalflare-rename-resname $SAVE_DIR/7-$PDB_ID-resnames.pdb TIP WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb
metalflare-rename-resname $SAVE_DIR/7-$PDB_ID-resnames.pdb TIP3 WAT --output $SAVE_DIR/7-$PDB_ID-resnames.pdb
```

## pdb4amber

```bash
pdb4amber -i $SAVE_DIR/7-$PDB_ID-resnames.pdb > $SAVE_DIR/8-$PDB_ID-pdb4amber.pdb
```

## Final

```bash
cp $SAVE_DIR/8-$PDB_ID-pdb4amber.pdb $SAVE_DIR/$PDB_ID-final.pdb
```

<!-- LINKS -->

[pdb2pqr]: https://github.com/Electrostatics/pdb2pqr
