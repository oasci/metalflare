# 001 - roGFP2 simulation

This experiment performs classical MD simulations of the [reduced form of roGFP2](../../methods/01-protocols/gfp-definitions.md#reduced-form) for benchmarking structural changes of [mseGFP mutations](../../methods/01-protocols/gfp-definitions.md#mseGFP).

## Protein preparation

We prepared the roGFP2 structure (PDB: `1JC0`) with the [protein preparation from RCSB protocol][protocol-protein-prep].
The [final structure][final-pdb] is shown below with labels on the experiment's relevant residues.

<div id="prepped-pdb-view" class="mol-container"></div>
<script>
var viewer1 = $3Dmol.createViewer(
    document.querySelector('#prepped-pdb-view'), { backgroundAlpha: '0.0' }
);
var pdbUri = '../../data/001-rogfp-md/structures/protein/1JC0-final.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer1.addModel( data, 'pdb' );
        viewer1.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer1.setStyle({resn: 'CRO'}, {stick: {}});
        viewer1.setStyle({resi: 145}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer1.setStyle({resi: 202}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer1.addLabel(
            "CRO 65",
            {screenOffset: new $3Dmol.Vector2(0, 0), backgroundOpacity: 0.8},
            {resi: 65}, false
        )
        viewer1.addLabel(
            "CYM 145",
            {screenOffset: new $3Dmol.Vector2(-100, 20), backgroundOpacity: 0.8},
            {resi: 145}, false
        )
        viewer1.addLabel(
            "CYM 202",
            {screenOffset: new $3Dmol.Vector2(30, 20), backgroundOpacity: 0.8},
            {resi: 202}, false
        )
        viewer1.zoomTo();
        viewer1.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

!!! warning

    Because the chromophore (`CRO 65`) is a non-standard residue, [PDB2PQR][protocol-protein-prep-pdb2pqr] does not place any hydrogens and are missing in [1JC0-final.pdb][final-pdb].
    Hydrogen atoms are added in the next protocol (tleap).

<!-- LINKS -->

[protocol-protein-prep]: ../../methods/02-protein-prep.md
[protocol-protein-prep-pdb2pqr]: ../../methods/02-protein-prep.md#protonation-and-steric-clashes
[final-pdb]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/data/001-rogfp-md/structures/protein/1JC0-final.pdb?ref_type=heads
