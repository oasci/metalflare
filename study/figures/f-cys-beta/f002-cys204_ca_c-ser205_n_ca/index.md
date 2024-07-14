# f002-cys204_ca_c-ser205_n_ca

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
        viewer.addLabel("CA", {}, {chain: "A", resi: 204, atom: "CA"})
        viewer.addLabel("C", {}, {chain: "A", resi: 204, atom: "C"})
        viewer.addLabel("N", {}, {chain: "A", resi: 205, atom: "N"})
        viewer.addLabel("CA", {}, {chain: "A", resi: 205, atom: "CA"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -182.16736190247866, -8.725446987364373, -41.92004156302197, 105.89065560538575, 0.1559475855759166, -0.9340336582546589, 0.16705031603383494, 0.2745098681531451 ]);
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
![](./f002-cys204_ca_c-ser205_n_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f002-cys204_ca_c-ser205_n_ca/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./f002-cys204_ca_c-ser205_n_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f002-cys204_ca_c-ser205_n_ca/pmf-info.md"
