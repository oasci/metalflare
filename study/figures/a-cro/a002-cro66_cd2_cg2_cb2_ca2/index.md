# a002-cro66_cd2_cg2_cb2_ca2

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
        viewer.setStyle({chain: 'A', resi: 146}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 147}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 148}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 204}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 205}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.addLabel("CD2", {}, {chain: "A", resi: resi1, atom: "CD2"})
        viewer.addLabel("CG2", {}, {chain: "A", resi: resi1, atom: "CG2"})
        viewer.addLabel("CB2", {}, {chain: "A", resi: resi1, atom: "CB2"})
        viewer.addLabel("CA2", {}, {chain: "A", resi: resi1, atom: "CA2"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -176.11118801308325, -9.101064812111433, -43.02812016330148, 105.68603275728341, 0.15291371871933526, -0.8780426305553863, -0.4533077704476863, -0.013061347807533916 ]);
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

## Probability density function

<figure markdown>
![](./a002-cro66_cd2_cg2_cb2_ca2-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a001-cro66_og1_cb1_ca1_c1/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./a002-cro66_cd2_cg2_cb2_ca2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a001-cro66_og1_cb1_ca1_c1/pmf-info.md"
