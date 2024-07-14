# a003-cro66_cg2_cb2_ca2_c2

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
        viewer.addLabel("CG2", {}, {chain: "A", resi: resi1, atom: "CG2"})
        viewer.addLabel("CB2", {}, {chain: "A", resi: resi1, atom: "CB2"})
        viewer.addLabel("CA2", {}, {chain: "A", resi: resi1, atom: "CA2"})
        viewer.addLabel("C2", {}, {chain: "A", resi: resi1, atom: "C2"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -176.32819422335245, -11.5435395679008, -43.70052957148548, 113.677239617803, 0.1981859208096388, -0.8621311737597761, -0.46627736602015957, 0.0061317176816809565 ]);
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
![](./a003-cro66_cg2_cb2_ca2_c2-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a003-cro66_cg2_cb2_ca2_c2/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./a003-cro66_cg2_cb2_ca2_c2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a003-cro66_cg2_cb2_ca2_c2/pmf-info.md"
