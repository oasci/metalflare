# B004: Cys147 CA - Cys204 CA

## Probability density function

<figure markdown>
![](./b004-cys147_ca-cys204_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cys/b004-cys147_ca-cys204_ca/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./b004-cys147_ca-cys204_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cys/b004-cys147_ca-cys204_ca/pmf-info.md"

## Visualization

<div id="b004-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('b004-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/f001.molj", "molj");
    });
});
</script>
