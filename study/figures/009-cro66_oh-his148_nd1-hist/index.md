# 009-cro66_oh-his148_nd1-hist

TODO:

<div id="rogfp-view" class="mol-container"></div>
<script>
var uri = 'https://files.rcsb.org/view/1jc0.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#rogfp-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 66;
        let atom1Name = "OH";
        let resi2 = 148;
        let atom2Name = "ND1";
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: resi1}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: resi2}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -182.846002868459, -11.360836858450071, -44.309189914886645, 111.40397984824831, -0.6066648911438127, 0.27532366621754684, -0.6413719126788533, 0.3805215608880305 ]);
        let atom1 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi1, atom: atom1Name}
        )[0];
        let atom2 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi2, atom: atom2Name}
        )[0];
        viewer.addCylinder(
            {
                dashed: true,
                start: {x: atom1.x, y: atom1.y, z: atom1.z},
                end: {x: atom2.x, y: atom2.y, z: atom2.z},
                radius: 0.1,
                color: "#00b4d8"
            }
        );
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
![](./009-cro66_oh-his148_nd1-hist.svg)
</figure>
