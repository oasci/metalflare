# GFP definitions

## eGFP

[Enhanced GFP][2y0g] (eGFP), first introduced by [Heim et al.][egfp paper], has `S65T` and `F64L` mutations from the wild type protein.

<div id="egfp-view" class="mol-container"></div>
<script>
var viewer1 = $3Dmol.createViewer(
    document.querySelector('#egfp-view'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/2y0g.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer1.addModel( data, 'pdb' );
        viewer1.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer1.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer1.setView([-24.438022254475037, -28.531666182873735, -38.325676826318364, -13.834375761483017, 0.28482918472442975, -0.7868184859871322, 0.21876734309052062, 0.5019261451999166]);
        viewer1.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

## roGFP

The [redox-sensitive GFP][rogfp paper] (roGFP) is derived from [eGFP](#egfp) with two additional mutations: `S147C` and `Q204C`.
Introduced as roGFP2 [Hanson et al.][rogfp paper], this forms a reversible formation of a [reduced][1jc0] and [oxidized][1jc1] disulfide bridge between `147` and `204`.

!!! note

    Due to renumbering residues in our [protein preparation pipeline](../02-protein-prep.md), these residues are `145` and `202` in our simulations.

### Reduced form

[1JC0][1jc0] shows the reduced (i.e., broken) form of `147`-`204` disulfide bond.

<div id="rogfp-reduced-view" class="mol-container"></div>
<script>
var viewer2 = $3Dmol.createViewer(
    document.querySelector('#rogfp-reduced-view'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/1jc0.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer2.addModel( data, 'pdb' );
        viewer2.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer2.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer2.setStyle({chain: 'A', resi: '147'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer2.setStyle({chain: 'A', resi: '204'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer2.setStyle({chain: 'B'}, {});
        viewer2.setStyle({chain: 'C'}, {});
        viewer2.setView([ -185.57030645037585, -7.386171514064613, -43.82421102060345, 88.52221071593, -0.2662199084146719, 0.4273867540214552, -0.7250566897044959, -0.4698513802953507 ]);
        viewer2.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

### Oxidized form

[1JC1][1jc1] shows the oxidized (i.e., formed) form of `147`-`204` disulfide bond.

<div id="rogfp-oxidized-view" class="mol-container"></div>
<script>
var viewer3 = $3Dmol.createViewer(
    document.querySelector('#rogfp-oxidized-view'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/1JC1.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer3.addModel( data, 'pdb' );
        viewer3.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer3.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer3.setStyle({chain: 'A', resi: '147'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer3.setStyle({chain: 'A', resi: '204'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer3.setStyle({chain: 'B'}, {});
        viewer3.setStyle({chain: 'C'}, {});
        viewer3.setView([ -186.81094603794497, -7.433604180840877, -44.03122729199806, 88.61698273769923, -0.28983518838466366, 0.46421243736754847, -0.6828026754229519, -0.4840277709001468 ]);
        viewer3.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

## mseGFP

[Metal-sensing GFP (mseGFP)][8dta] is similar to roGFP2, but mutates `S147C` and `S202C` from [eGFP][2y0g].
Since roGFP2 can no longer make the `147`-`204` disulfide bond, all mseGFP simulations will start from the [reduced form](#reduced-form) with `C204S` and `S202C` mutations.

!!! note

    Due to renumbering residues in our [protein preparation pipeline](../02-protein-prep.md), these residues are `145` and `200` in our simulations.

<div id="msegfp-view" class="mol-container"></div>
<script>
var viewer4 = $3Dmol.createViewer(
    document.querySelector('#msegfp-view'), { backgroundAlpha: '0.0' }
);
var pdbUri = 'https://files.rcsb.org/view/8DTA.pdb';
jQuery.ajax( pdbUri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        viewer4.addModel( data, 'pdb' );
        viewer4.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer4.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer4.setStyle({chain: 'A', resi: '147'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer4.setStyle({chain: 'A', resi: '202'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer4.setView([ -60.64682338153259, -20.114962159611807, 0.5702077286702113, 80.5194132281471, -0.15077826938374425, 0.19679882644092048, -0.8102144809849335, -0.5311201654949984 ]);
        viewer4.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
// function getState() {
//   console.log(console.log(viewer4.getView()));
// };
</script>
<!-- <button onclick="getState()">Click me to get 3DMol.js state</button> -->

<!-- LINKS -->

[egfp paper]: https://doi.org/10.1038/373663b0
[1jc0]: https://www.rcsb.org/structure/1jc0
[1jc1]: https://www.rcsb.org/structure/1jc1
[rogfp paper]: https://doi.org/10.1074/jbc.M312846200
[2y0g]: https://www.rcsb.org/structure/2y0g
[8dta]: https://www.rcsb.org/structure/8DTA
