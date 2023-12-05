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
        viewer1.zoomTo({chain: 'A'});
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
        viewer2.zoomTo({chain: 'A', resi: '204'});
        viewer2.zoom(0.6);
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
        viewer3.zoomTo({chain: 'A', resi: '204'});
        viewer3.zoom(0.6);
        viewer3.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

## mseGFP

Metal-sensing GFP (mseGFP) is similar to roGFP2, but mutates `S147C` and `S202C` from [eGFP][2y0g].
Since roGFP2 can no longer make the `147`-`204` disulfide bond, all mseGFP simulations will start from the [reduced form](#reduced-form) with `C202S` and `S202C` mutations.

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
        viewer4.zoomTo({chain: 'A', resi: '202'});
        viewer4.zoom(0.6);
        viewer4.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
});
</script>

<!-- LINKS -->

[egfp paper]: https://doi.org/10.1038/373663b0
[1jc0]: https://www.rcsb.org/structure/1jc0
[1jc1]: https://www.rcsb.org/structure/1jc1
[rogfp paper]: https://doi.org/10.1074/jbc.M312846200
[2y0g]: https://www.rcsb.org/structure/2y0g
