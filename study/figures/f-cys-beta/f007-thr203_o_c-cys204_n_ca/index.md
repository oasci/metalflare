# f007-thr203_o_c-cys204_n_ca

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
        let resi1 = 204;
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
        viewer.addLabel("O", {}, {chain: "A", resi: 203, atom: "O"})
        viewer.addLabel("C", {}, {chain: "A", resi: 203, atom: "C"})
        viewer.addLabel("N", {}, {chain: "A", resi: 204, atom: "N"})
        viewer.addLabel("CA", {}, {chain: "A", resi: 204, atom: "CA"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -183.67039134673692, -8.790563885920129, -42.01344809829062, 111.72365444332763, 0.0934415604548533, -0.24428734671013674, 0.8726447904804752, 0.41241173196742426 ]);
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
![](./f007-thr203_o_c-cys204_n_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f007-thr203_o_c-cys204_n_ca/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./f007-thr203_o_c-cys204_n_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f007-thr203_o_c-cys204_n_ca/pmf-info.md"
