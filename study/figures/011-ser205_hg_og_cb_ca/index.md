# 011-ser205_hg_og_cb_ca

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
        viewer.setStyle({chain: 'A', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 204}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 205}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.addLabel("HG", {screenOffset: {x: -20, y: 50}}, {chain: "A", resi: resi1, atom: "OG"})
        viewer.addLabel("OG", {}, {chain: "A", resi: resi1, atom: "OG"})
        viewer.addLabel("CB", {}, {chain: "A", resi: resi1, atom: "CB"})
        viewer.addLabel("CA", {}, {chain: "A", resi: resi1, atom: "CA"})
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -178.82870089204013, -9.020158654661756, -44.60591775213842, 97.41737754014196, -0.2158946699173421, 0.9348675629168328, 0.14514172064286476, -0.24154919216600987 ]);
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
![](./011-ser205_hg_og_cb_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/011-ser205_hg_og_cb_ca/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./011-ser205_hg_og_cb_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/011-ser205_hg_og_cb_ca/pmf-info.md"
