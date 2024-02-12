# 013-thr203_o_c_ca_cb

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
        let resi1 = 203;
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
        viewer.addLabel("O", {}, {chain: "A", resi: resi1, atom: "O"})
        viewer.addLabel("C", {}, {chain: "A", resi: resi1, atom: "C"})
        viewer.addLabel("CA", {}, {chain: "A", resi: resi1, atom: "CA"})
        viewer.addLabel("CB", {}, {chain: "A", resi: resi1, atom: "CB"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -181.1940983641914, -10.88805946603958, -42.57839978149105, 105.68603275728339, 0.4498289015763534, 0.0013940082009127504, -0.11073459989002167, 0.8862222432521896 ]);
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
![](./013-thr203_o_c_ca_cb-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/013-thr203_o_c_ca_cb/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./013-thr203_o_c_ca_cb-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/013-thr203_o_c_ca_cb/pmf-info.md"
