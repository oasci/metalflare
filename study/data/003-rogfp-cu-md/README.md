# 003 - roGFP2 with Cu(I) simulation

This experiment performs classical MD simulations of the [reduced form of roGFP2](../../methods/01-protocols/gfp-definitions.md#reduced-form) with Cu(I) bound to two reduced cysteines.

## Placing Cu(I)

<div id="cu-placed-view" class="mol-container"></div>
<script>
var uri = './structures/protein/1JC0-Cu.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#cu-placed-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModelsAsFrames(data, "pdb");
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer.setStyle(
            {chain: 'A', resi: '145'}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.95}}
        );
        viewer.setStyle(
            {chain: 'A', resi: '146'}, {cartoon: {color: "spectrum", opacity: 0.95}}
        );
        viewer.setStyle(
            {chain: 'A', resi: '202'}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.95}}
        );
        viewer.setStyle(
            {chain: 'A', resi: '203'}, {cartoon: {color: "spectrum", opacity: 0.95}}
        );
        viewer.setStyle({resn: 'CU1'}, {sphere: {radius: 0.9, color: "0xa52a2a"}});
        viewer.setView([ -36.9147351738292, -38.79385525105598, -34.28640214794895, 72.28733464747603, 0.2666368765879281, -0.19611494476442937, 0.7836729219496249, 0.5256428976847463 ]);
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
        viewer.setStyle({resn: 'CU1'}, {sphere: {radius: 0.9, color: "0xa52a2a"}});
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
