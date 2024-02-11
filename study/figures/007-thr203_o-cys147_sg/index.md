# 007-thr203_o-cys147_sg

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
        let atom1Name = "O";
        let resi2 = 147;
        let atom2Name = "SG";
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
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -183.8551549370776, -8.180688339796237, -44.32967895720527, 102.98523157925786, 0.04006695437297919, -0.9920895090847902, 0.11566121531051957, 0.027848310614215 ]);
        let atom1 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi1, atom: atom1Name}
        )[0];
        let atom2 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi2, atom: atom2Name}
        )[0];
        viewer.addCylinder(
            {
                dashed: true,
                start: {x: atom1.x, y: atom1.y, z: atom1.z},
                end: {x: atom2.x, y: atom2.y, z: atom2.z},
                radius: 0.1,
                color: "#00b4d8"
            }
        );
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
![](./007-thr203_o-cys147_sg-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/007-thr203_o-cys147_sg/pdf-info.md"

### Bandwidth validation

<figure markdown>
![](./007-thr203_o-cys147_sg-hist.svg)
</figure>


## Potential of mean force

<figure markdown>
![](./007-thr203_o-cys147_sg-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/007-thr203_o-cys147_sg/pmf-info.md"
