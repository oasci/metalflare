# a001-cro66_og1_cb1_ca1_c1

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
        viewer.addLabel("OG1", {}, {chain: "A", resi: resi1, atom: "OG1"})
        viewer.addLabel("CB1", {}, {chain: "A", resi: resi1, atom: "CB1"})
        viewer.addLabel("CA1", {}, {chain: "A", resi: resi1, atom: "CA1"})
        viewer.addLabel("C1", {}, {chain: "A", resi: resi1, atom: "C1"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -178.54845283019642, -8.287854743737649, -41.32961916913211, 96.76437948169743, 0.12322892294046181, -0.986881068950569, 0.04557929448433187, 0.0938238573723442 ]);
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
![](./a001-cro66_og1_cb1_ca1_c1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a001-cro66_og1_cb1_ca1_c1/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./a001-cro66_og1_cb1_ca1_c1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a001-cro66_og1_cb1_ca1_c1/pmf-info.md"
