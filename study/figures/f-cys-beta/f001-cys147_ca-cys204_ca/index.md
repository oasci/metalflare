# f001-cys147_ca-cys204_ca

TODO:

<div id="f001-view" class="mol-container"></div>
<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('f001-view', {
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
        // viewer.loadStructureFromUrl("/analysis/005-rogfp-glh-md/data/traj/frame_7722.pdb", "pdb");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/f001.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./f001-cys147_ca-cys204_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f001-cys147_ca-cys204_ca/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./f001-cys147_ca-cys204_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f001-cys147_ca-cys204_ca/pmf-info.md"