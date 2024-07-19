# f004-ser205_h-asn146_o

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
        let atom1Name = "N";
        let resi2 = 146;
        let atom2Name = "O";
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
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -180.33196291381824, -12.112329257963113, -46.47604478199624, 97.49843708210796, 0.17474428971354652, -0.9161577310617698, 0.2631075810734197, 0.2467667842517564 ]);
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
![](./f004-ser205_h-asn146_o-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f004-ser205_h-asn146_o/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./f004-ser205_h-asn146_o-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f004-ser205_h-asn146_o/pmf-info.md"
