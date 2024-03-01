# 008-thr203_o-cys204_sg

TODO:

<div id="rogfp-view" class="mol-container"></div>
<script>
var uri = '../../analysis/987-rogfp-cu-md-select/pdbs/single-select/single-select-3-rogfp_cu.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#rogfp-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 201;
        let atom1Name = "O";
        let resi2 = 202;
        let atom2Name = "SG";
        viewer.addModel( data, 'pdb' );
        viewer.addLabel(
            "3.15 Å",
            {screenOffset: new $3Dmol.Vector2(-5, 150), backgroundOpacity: 0.8},
            {resi: resi1}, false
        )
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
        viewer.setView([ -34.98128057494662, -53.51902927276367, -49.26445180537328, 112.63787023958388, 0.5442619345726372, -0.3702171717399282, 0.7235920071257489, 0.20768437482841262 ]);
        let atom1 = viewer.getModel().selectedAtoms(
            {chain: 'X', resi: resi1, atom: atom1Name}
        )[0];
        let atom2 = viewer.getModel().selectedAtoms(
            {chain: 'X', resi: resi2, atom: atom2Name}
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

<div id="rogfp-view-2" class="mol-container"></div>
<script>
var uri = '../../analysis/987-rogfp-cu-md-select/pdbs/single-select/single-select-2-rogfp.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#rogfp-view-2'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 201;
        let atom1Name = "O";
        let resi2 = 202;
        let atom2Name = "SG";
        viewer.addModel( data, 'pdb' );
        viewer.addLabel(
            "3.56 Å",
            {screenOffset: new $3Dmol.Vector2(40, 110), backgroundOpacity: 0.8},
            {resi: resi1}, false
        )
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
        viewer.setView( [ -46.894827097138084, -28.578266761678005, -35.329224280794506, 103.78068554565341, 0.29973746457542993, -0.135574783397343, 0.9046901717885979, 0.27076303940119467 ]);
        let atom1 = viewer.getModel().selectedAtoms(
            {chain: 'X', resi: resi1, atom: atom1Name}
        )[0];
        let atom2 = viewer.getModel().selectedAtoms(
            {chain: 'X', resi: resi2, atom: atom2Name}
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
![](./008-thr203_o-cys204_sg-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/008-thr203_o-cys204_sg/pdf-info.md"

### Bandwidth validation

<figure markdown>
![](./008-thr203_o-cys204_sg-hist.svg)
</figure>

## Potential of mean force

<figure markdown>
![](./008-thr203_o-cys204_sg-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/008-thr203_o-cys204_sg/pmf-info.md"
