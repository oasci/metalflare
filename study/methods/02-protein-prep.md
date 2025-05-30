# 02 - Protein preparation from RCSB

Crystallographic structures from the Protein Data Bank require extensive preprocessing before use in molecular dynamics simulations. Raw PDB files contain artifacts, missing atoms, and structural inconsistencies that must be addressed systematically. Our pipeline consists of five distinct stages, each handling specific aspects of structure preparation and simulation setup.

The complete workflow proceeds through: (1) protein structure cleaning and optimization, (2) force field parameterization and solvation, (3) energy minimization preparation, (4) equilibration protocol setup, and (5) production simulation configuration. Each stage builds upon the previous one, ensuring reproducibility and maintaining careful control over simulation conditions.

### Initial Structure Acquisition and Processing

The preparation begins with downloading the target structure (1JC0) from the RCSB Protein Data Bank. This particular structure represents the reduced form of roGFP2, which serves as our reference state for Cu(I) binding investigations. The crystallographic structure contains multiple chains, solvent molecules, and crystallization artifacts that require systematic removal or modification.

```bash
export PDB_ID="1JC0"
wget https://files.rcsb.org/download/$PDB_ID.pdb -O structures/protein/0-$PDB_ID.pdb
```

The initial download yields the complete crystallographic asymmetric unit, including multiple protein chains and extensive crystallographic waters.


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


For our simulation purposes, we focus exclusively on chain A, which contains the complete roGFP2 structure with the characteristic chromophore and cysteine residues critical for Cu(I) coordination.

### Chain Selection and Water Retention

Chain selection removes extraneous protein chains while preserving crystallographic waters that may play important structural roles. The `metalflare-select-atoms` utility employs MDAnalysis selection syntax to isolate chain A while maintaining all solvent molecules:

```bash
metalflare-select-atoms 0-1JC0.pdb 1-1JC0-chain-A.pdb --select_str "chainID A"
```

We will only be working with one protein and need to select chain `A`.
We use the `metalflare-select-atoms` script which just drives `pdb.select.run_select_atoms()`.

```bash
metalflare-select-atoms $SAVE_DIR/0-$PDB_ID.pdb $SAVE_DIR/1-$PDB_ID-chain-A.pdb --select_str "chainID A"
```

We keep all crystallographic waters so the CRO-coordinating water's position is maintained.

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
        viewer.setStyle({resn: 'HOH'}, {sphere: {scale: '0.3', opacity: '0.95'}});
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

This selective approach is crucial because crystallographic waters often occupy functionally important positions, particularly those coordinating with the chromophore or stabilizing key structural elements. Removing these waters indiscriminately can create artificial cavities or destabilize important hydrogen bonding networks during subsequent simulation.

### Structure Filtering and Coordinate Optimization

The filtering step removes non-essential PDB records, retaining only `ATOM` and `HETATM` entries necessary for molecular dynamics. This cleanup eliminates crystallographic metadata, experimental annotations, and other records that can interfere with downstream force field assignment:

```bash
metalflare-filter-pdb 1-1JC0-chain-A.pdb --output 2-1JC0-filtered.pdb
```

Following filtration, the structure undergoes geometric centering to position the protein at the coordinate system origin. This centering simplifies subsequent box construction and ensures consistent placement across different simulation systems:

```bash
metalflare-center 2-1JC0-filtered.pdb --output 3-1JC0-centered.pdb
```

### Rotational Optimization for Solvation Efficiency

The rotational optimization step represents a critical but often overlooked aspect of simulation preparation. The `metalflare-minimize-box` utility systematically tests different protein orientations to identify the configuration that minimizes the required solvation volume:

```bash
metalflare-minimize-box 3-1JC0-centered.pdb --output 4-1JC0-rotated.pdb
```

This reduction significantly decreases computational requirements while maintaining proper solvation around all protein surfaces. The algorithm evaluates multiple rotational configurations and selects the orientation that achieves the most compact solvation shell without creating buried surfaces or artificial contacts.

### Residue Numbering and Chemical State Assignment

Crystallographic structures often contain inconsistent residue numbering or non-standard residue designations that must be standardized for force field compatibility. The unification process ensures sequential residue numbering and consistent chain identifiers:

```bash
metalflare-unify-resids 4-1JC0-rotated.pdb --output 5-1JC0-residues.pdb
```

### Chemical Modifications and Protonation State Management

#### Selenomethionine Conversion

Crystallographic structures frequently contain selenomethionine (MSE) residues introduced during crystallization for phasing purposes. These artificial modifications must be reverted to native methionine residues for accurate simulation:

```bash
metalflare-rename-resname 5-1JC0-residues.pdb MSE MET --output 5-1JC0-residues.pdb
```

#### Cysteine State Assignment

The roGFP2 structure contains two critical cysteine residues (Cys147 and Cys204) that coordinate Cu(I) binding. These residues must be properly configured in their reduced, thiol-bearing state to enable metal coordination. The deprotonated cysteine state (CYM) is assigned to these specific residues:

```bash
metalflare-rename-resname 5-1JC0-residues.pdb CYS CYM --include 147 204 --output 5-1JC0-residues.pdb
```

This assignment is crucial because standard cysteine residues (CYS) model the protonated thiol form, which cannot coordinate metals effectively. The deprotonated cysteine state (CYM) provides the anionic sulfur atoms necessary for Cu(I) coordination chemistry.

#### Glutamic Acid Protonation

Specific glutamic acid residues may require protonation state adjustment based on local electrostatic environments and pH considerations. Residue 220 is converted to the protonated form (GLH) to maintain proper charge balance and local structural stability:

```bash
metalflare-rename-resname 5-1JC0-residues.pdb GLU GLH --include 220 --output 5-1JC0-residues.pdb
```

### Automated Protonation State Prediction

The PDB2PQR program provides automated protonation state assignment for ionizable residues based on local electrostatic environments and pH conditions. This tool addresses the fundamental challenge that crystallographic structures lack hydrogen atoms, which are essential for accurate force field modeling:

```bash
pdb2pqr --log-level INFO --ff=AMBER --keep-chain --ffout=AMBER 5-1JC0-residues.pdb 6-1JC0-pdb2pqr.pdb
```

PDB2PQR employs sophisticated algorithms to predict optimal protonation states for histidine, aspartic acid, glutamic acid, and lysine residues. The program considers local electrostatic interactions, hydrogen bonding opportunities, and bulk pH conditions to determine the most probable protonation pattern. However, PDB2PQR cannot process non-standard residues like the GFP chromophore, requiring manual intervention for these special cases.

The merge operation combines PDB2PQR output with the original structure to recover any atoms that PDB2PQR could not process:

```bash
metalflare-merge-pdbs 6-1JC0-pdb2pqr.pdb 5-1JC0-residues.pdb --output 6-1JC0-pdb2pqr.pdb
```

### Water Molecule Standardization

Crystallographic structures may contain water molecules with various naming conventions (HOH, TIP, TIP3). Standardization to WAT ensures consistent recognition by Amber force field tools:

```bash
metalflare-rename-resname 6-1JC0-pdb2pqr.pdb HOH WAT --output 7-1JC0-resnames.pdb
metalflare-rename-resname 7-1JC0-resnames.pdb TIP WAT --output 7-1JC0-resnames.pdb
metalflare-rename-resname 7-1JC0-resnames.pdb TIP3 WAT --output 7-1JC0-resnames.pdb
metalflare-unify-waters 7-1JC0-resnames.pdb --output 7-1JC0-resnames.pdb
```

The unification process also ensures that water molecules have consistent atom naming and connectivity, preventing force field assignment errors during subsequent parameterization steps.

### Final Amber Compatibility Processing

The pdb4amber utility performs final compatibility checks and modifications to ensure full compatibility with Amber force field tools:

```bash
pdb4amber -i 7-1JC0-resnames.pdb > 8-1JC0-pdb4amber.pdb 2> pdb4amber.err
```

This tool addresses Amber-specific requirements such as terminal residue modifications, special residue recognition, and atom naming conventions that differ from standard PDB formatting. The processed structure becomes the final input for force field parameterization and system construction.

<!-- LINKS -->

[pdb2pqr]: https://github.com/Electrostatics/pdb2pqr
