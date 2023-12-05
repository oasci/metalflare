# 002 - Cu(I) positioning

TODO:

## Subsystem

TODO:

<div id="subsystem-initial-view" class="mol-container"></div>
<script>
var uri = './structures/01-initial/subsystem.xyz';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#subsystem-initial-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'xyz' );
        viewer.setStyle({}, {stick: {}});
        viewer.setView([ -41.294419203821676, -35.5929660191083, -35.54826140127385, 87.78780484222388, -0.23736539886223026, 0.34657018820359475, -0.7530590848177492, -0.5064077278685156 ]);
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + uri + ": " + err );
    },
});
</script>

## Solvent optimization

<div id="geo-opt-solv-view" class="mol-container"></div>
<script>
var uri = './calculations/01-solvent-opt/xtbopt.log';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#geo-opt-solv-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModelsAsFrames(data, "xyz");
        viewer.animate({interval: 150, loop: "forward", reps: 0});
        viewer.setStyle({}, {stick: {}});
        viewer.setView([ -41.294419203821676, -35.5929660191083, -35.54826140127385, 87.78780484222388, -0.23736539886223026, 0.34657018820359475, -0.7530590848177492, -0.5064077278685156 ]);
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

## Cu(I) placement

<div id="copper-placement-view" class="mol-container"></div>
<script>
var uri = './calculations/02-dock-copper/xtbopt.log';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#copper-placement-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModelsAsFrames(data, "xyz");
        viewer.animate({interval: 150, loop: "forward", reps: 0});
        viewer.setStyle({}, {stick: {}});
        viewer.setStyle({elem: "Cu"}, {sphere: {scale: 0.5}});
        viewer.setView([ -47.37641165022477, -34.81989291611326, -34.415538199484196, 107.15776591487352, -0.17350995400483246, 0.18953312050907078, -0.7973543906475701, -0.5460745991221131 ]);
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>
