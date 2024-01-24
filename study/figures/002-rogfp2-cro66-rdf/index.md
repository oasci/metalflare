# 002-rogfp2-cro66-rdf

TODO:

## O rdf

TODO:

<div id="goo-view" class="mol-container"></div>
<script>
var uri = '../../../data/003-rogfp-cu-md/structures/protein/1JC0-Cu.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#goo-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 65;
        let atom1Name = "OH";
        let dist = 8;
        viewer.addModelsAsFrames(data, "pdb");
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.6}});
        viewer.setStyle(
            {chain: 'A', elem: "O", within: {distance: dist, sel: {chain: 'A', resi: resi1, atom: atom1Name}}},
            {sphere: {radius: 0.3}}
        );
        viewer.setStyle({chain: 'A', resi: resi1}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setView([ -40.27769001496433, -40.93068228479875, -39.21748927182831, 91.73385345736014, -0.8716373326314918, -0.2627468253148699, 0.12873747533973695, -0.393241819486643 ]);
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


<figure markdown>
![](./cro66-rdf-o.svg)
</figure>

## N rdf

<figure markdown>
![](./cro66-rdf-n.svg)
</figure>

<div id="gon-view" class="mol-container"></div>
<script>
var uri = '../../../data/003-rogfp-cu-md/structures/protein/1JC0-Cu.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#gon-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 65;
        let atom1Name = "OH";
        let dist = 8;
        viewer.addModelsAsFrames(data, "pdb");
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.6}});
        viewer.setStyle(
            {chain: 'A', elem: "N", within: {distance: dist, sel: {chain: 'A', resi: resi1, atom: atom1Name}}},
            {sphere: {radius: 0.3}}
        );
        viewer.setStyle({chain: 'A', resi: resi1}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setView([ -40.27769001496433, -40.93068228479875, -39.21748927182831, 91.73385345736014, -0.8716373326314918, -0.2627468253148699, 0.12873747533973695, -0.393241819486643 ]);
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

## H rdf

<figure markdown>
![](./cro66-rdf-h.svg)
</figure>
