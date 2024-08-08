# B005: Cys147 SG - Cys204 SG

## Probability density function

<figure markdown>
![](./b005-cys147_sg-cys204_sg-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cys/b005-cys147_sg-cys204_sg/pdf-info.md"

## Visualization

<div id="b005-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('b005-view', {
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
