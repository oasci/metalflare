# GFP definitions

## eGFP

[Enhanced GFP][2y0g] (eGFP), first introduced by [Heim et al.][egfp paper], has `S65T` and `F64L` mutations from the wild type protein.

!!! quote "2Y0G"
    <div id="2Y0G-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('2Y0G-view', {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: true,
        viewportShowAnimation: false,
        pdbProvider: 'rcsb',
    }).then(viewer => {
        // viewer.loadPdb("2Y0G");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/2Y0G.molj", "molj");
    });
});
</script>

## roGFP

The [redox-sensitive GFP][rogfp paper] (roGFP) is derived from [eGFP](#egfp) with two additional mutations: `S147C` and `Q204C`.
Introduced as roGFP2 [Hanson et al.][rogfp paper], this forms a reversible formation of a [reduced][1jc0] and [oxidized][1jc1] disulfide bridge between `147` and `204`.

!!! note

    Due to renumbering residues in our [protein preparation pipeline](../02-protein-prep.md), these residues are `145` and `202` in our simulations.

### Reduced form

[1JC0][1jc0] shows the reduced (i.e., broken) form of `147`-`204` disulfide bond.

!!! quote "1JC0"
    <div id="1JC0-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('1JC0-view', {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: true,
        viewportShowAnimation: false,
        pdbProvider: 'rcsb',
    }).then(viewer => {
        // viewer.loadPdb("1JC0");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/1JC0.molj", "molj");
    });
});
</script>

### Oxidized form

[1JC1][1jc1] shows the oxidized (i.e., formed) form of `147`-`204` disulfide bond.

!!! quote "1JC1"
    <div id="1JC1-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('1JC1-view', {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: true,
        viewportShowAnimation: false,
        pdbProvider: 'rcsb',
    }).then(viewer => {
        // viewer.loadPdb("1JC1");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/1JC1.molj", "molj");
    });
});
</script>

## mseGFP

[Metal-sensing GFP (mseGFP)][8dta] is similar to roGFP2, but mutates `S147C` and `S202C` from [eGFP][2y0g].
Since roGFP2 can no longer make the `147`-`204` disulfide bond, all mseGFP simulations will start from the [reduced form](#reduced-form) with `C204S` and `S202C` mutations.

!!! note

    Due to renumbering residues in our [protein preparation pipeline](../02-protein-prep.md), these residues are `145` and `200` in our simulations.
<!--
<div id="msegfp-view" class="mol-container"></div>
<script>
var uri = 'https://files.rcsb.org/view/8DTA.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#msegfp-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({chain: 'A', resn: 'CRO'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer.setStyle({chain: 'A', resi: '147'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer.setStyle({chain: 'A', resi: '202'}, {stick: {}, cartoon: {color: "spectrum"}});
        viewer.setView([ -60.64682338153259, -20.114962159611807, 0.5702077286702113, 80.5194132281471, -0.15077826938374425, 0.19679882644092048, -0.8102144809849335, -0.5311201654949984 ]);
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script> -->

<!-- LINKS -->

[egfp paper]: https://doi.org/10.1038/373663b0
[1jc0]: https://www.rcsb.org/structure/1jc0
[1jc1]: https://www.rcsb.org/structure/1jc1
[rogfp paper]: https://doi.org/10.1074/jbc.M312846200
[2y0g]: https://www.rcsb.org/structure/2y0g
[8dta]: https://www.rcsb.org/structure/8DTA
