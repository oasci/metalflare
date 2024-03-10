# 004 - roGFP2-oxd simulation

This experiment performs classical MD simulations of the [oxidized form of roGFP2](../../methods/01-protocols/gfp-definitions.md#oxidized-form) for benchmarking structural changes of [mseGFP mutations](../../methods/01-protocols/gfp-definitions.md#msegfp).

## Protein preparation

We prepared the roGFP2 structure (PDB: `1JC1`) with the [protein preparation from RCSB protocol][protocol-protein-prep].
The [final structure][final-pdb] is shown below with labels on the experiment's relevant residues.

<div id="prepped-pdb-view" class="mol-container"></div>
<script>
var uri = '../../data/001-rogfp-md/structures/protein/1JC1-final.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#prepped-pdb-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({resn: 'CRO'}, {stick: {}});
        viewer.setStyle({resi: 145}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer.setStyle({resi: 202}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer.addLabel(
            "CRO 65",
            {screenOffset: new $3Dmol.Vector2(0, 0), backgroundOpacity: 0.8},
            {resi: 65}, false
        )
        viewer.addLabel(
            "CYM 145",
            {screenOffset: new $3Dmol.Vector2(-100, 20), backgroundOpacity: 0.8},
            {resi: 145}, false
        )
        viewer.addLabel(
            "CYM 202",
            {screenOffset: new $3Dmol.Vector2(30, 20), backgroundOpacity: 0.8},
            {resi: 202}, false
        )
        // viewer.setView([ -0.7561101750598701, -0.9271423446320399, 2.965827751298417, 49.265373924881985, 0.37232883239820697, -0.4757222855340383, 0.6628384467092744, 0.4423852858937514 ]);
        viewer.zoomTo({chain: "A"})
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

!!! warning

    Because the chromophore (`CRO 65`) is a non-standard residue, PDB2PQR does not place any hydrogens and are missing in [1JC0-final.pdb][final-pdb].
    Hydrogen atoms are added in the next protocol (tleap).

## Production

Each frame represents a stride of 2.5 ns for a total of 100 ns.

<div id="prod-npt-view" class="mol-container"></div>
<script>
var uri = './simulations/05-prod/run-01/outputs/08_prod_npt.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#prod-npt-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModelsAsFrames(data, "pdb");
        viewer.animate({interval: 200, loop: "forward", reps: 0});
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({resn: 'CRO'}, {stick: {}});
        viewer.setStyle({resi: 145}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer.setStyle({resi: 202}, {stick: {}, cartoon: {color: 'spectrum'}});
        viewer.setView([ -31.023800442233295, -32.70469651741289, -32.66362686567166, 9.517105134040186, -0.5129826573726642, 0.6134893978796151, -0.32219140994260353, -0.5066283127534278 ]);
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

<!-- LINKS -->

[protocol-protein-prep]: ../../methods/02-protein-prep.md
[final-pdb]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/data/004-rogfp-oxd/structures/protein/1JC1-final.pdb
