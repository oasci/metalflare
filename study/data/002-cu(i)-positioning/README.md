# 002 - Cu(I) positioning

TODO:

## Subsystem

TODO:

<div id="subsystem-initial-view" class="mol-container"></div>
<script>
var viewer1 = $3Dmol.createViewer(
    document.querySelector('#subsystem-initial-view'), { backgroundAlpha: '0.0' }
);
var uri = './structures/01-initial/subsystem.xyz';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer1.addModel( data, 'xyz' );
        viewer1.setStyle({}, {stick: {}});
        viewer1.setView([ -41.294419203821676, -35.5929660191083, -35.54826140127385, 87.78780484222388, -0.23736539886223026, 0.34657018820359475, -0.7530590848177492, -0.5064077278685156 ]);
        viewer1.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + uri + ": " + err );
    },
});
</script>

## Solvent optimization

<div id="geo-opt-solv-view" class="mol-container"></div>
<script>
var viewer2 = $3Dmol.createViewer(
    document.querySelector('#geo-opt-solv-view'), { backgroundAlpha: '0.0' }
);
var uri = './calculations/01-solvent-opt/xtbopt.log';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer2.addModelsAsFrames(data, "xyz");
        viewer2.animate({interval: 150, loop: "forward", reps: 0});
        viewer2.setStyle({}, {stick: {}});
        viewer2.setView([ -41.294419203821676, -35.5929660191083, -35.54826140127385, 87.78780484222388, -0.23736539886223026, 0.34657018820359475, -0.7530590848177492, -0.5064077278685156 ]);
        viewer2.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
function getState() {
  console.log(console.log(viewer2.getView()));
};
</script>
<!-- <button onclick="getState()">Click me to get 3DMol.js state</button> -->
