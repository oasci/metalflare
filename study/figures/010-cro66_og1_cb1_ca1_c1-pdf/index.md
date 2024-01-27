# 010-cro66_og1_cb1_ca1_c1-pdf

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
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 66}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 145}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 147}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 148}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 204}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 205}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.addLabel("OG1", {}, {chain: "A", resi: resi1, atom: "OG1"})
        viewer.addLabel("CB1", {}, {chain: "A", resi: resi1, atom: "CB1"})
        viewer.addLabel("CA1", {}, {chain: "A", resi: resi1, atom: "CA1"})
        viewer.addLabel("C1", {}, {chain: "A", resi: resi1, atom: "C1"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -177.60999010701107, -7.076870579532307, -41.60640597287334, 113.95634728539571, -0.20818863454283987, 0.9500304295785287, 0.0064224141737788815, -0.2325046836478871 ]);
        //viewer.zoomTo({chain: "A"})
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
![](./cro66_og1_cb1_ca1_c1-pdf.svg)
</figure>
