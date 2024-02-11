# 014-ser205_og_cb_ca_n

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
        let resi1 = 205;
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 66}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 145}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 147}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 148}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 204}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 205}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.addLabel("OG", {}, {chain: "A", resi: resi1, atom: "OG"})
        viewer.addLabel("CB", {}, {chain: "A", resi: resi1, atom: "CB"})
        viewer.addLabel("CA", {}, {chain: "A", resi: resi1, atom: "CA"})
        viewer.addLabel("N", {}, {chain: "A", resi: resi1, atom: "N"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -185.70717595348242, -9.400820059425605, -43.48492795815956, 113.90069785674157, 0.1345710483371672, -0.07152813837129907, 0.8996454582295998, 0.4091606137660668 ]);
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
![](./014-ser205_og_cb_ca_n-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/014-ser205_og_cb_ca_n/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./014-ser205_og_cb_ca_n-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/014-ser205_og_cb_ca_n/pmf-info.md"
