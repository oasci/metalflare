# f007-cys147_sg-cys204_sg

TODO:

<div id="f007-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('f007-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/f007.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./f007-cys147_sg-cys204_sg-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f007-cys147_sg-cys204_sg/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./f007-cys147_sg-cys204_sg-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/f007-cys147_sg-cys204_sg/pmf-info.md"