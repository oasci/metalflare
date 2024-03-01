# 012-thr203_hg1_og1_cb_cg2

TODO:

<div id="rogfp-view" class="mol-container"></div>
<script>
var uri = '../../analysis/987-rogfp-cu-md-select/pdbs/single-select/single-select-5-rogfp_cu.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#rogfp-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 201;
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'X'}, {cartoon: {color: 'spectrum', opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 65}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 143}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 144}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 145}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 146}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 201}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 202}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'X', resn: "CU1"}, {sphere: {radius: 1.0}});
        viewer.addLabel("HG1", {screenOffset: {x: -50, y: 50}}, {chain: "X", resi: resi1, atom: "OG1"})
        viewer.addLabel("OG1", {}, {chain: "X", resi: resi1, atom: "OG1"})
        viewer.addLabel("CB", {}, {chain: "X", resi: resi1, atom: "CB"})
        viewer.addLabel("CG2", {}, {chain: "X", resi: resi1, atom: "CG2"})
        viewer.setView([ -33.002337119885084, -34.727401945170485, -46.30887708298572, 113.35767842091929, -0.02381025599965196, 0.14284795244407691, -0.9838127028071448, 0.10554667207954437 ]);
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
![](./012-thr203_hg1_og1_cb_cg2-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/012-thr203_hg1_og1_cb_cg2/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./012-thr203_hg1_og1_cb_cg2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/012-thr203_hg1_og1_cb_cg2/pmf-info.md"
