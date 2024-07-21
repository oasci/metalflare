# l003-his148_hd1-asn146_o

TODO:

<div id="l003-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('l003-view', {
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
        // viewer.loadStructureFromUrl("/analysis/005-rogfp-glh-md/data/traj/frame_106403.pdb", "pdb");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/l003.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./l003-his148_hd1-asn146_o-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/l003-his148_hd1-asn146_o/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./l003-his148_hd1-asn146_o-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/l003-his148_hd1-asn146_o/pmf-info.md"
