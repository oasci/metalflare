# B006: Cys204_cb_ca_c_o

## Probability density function

<figure markdown>
![](./b006-cys204_cb_ca_c_o-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cys/b006-cys204_cb_ca_c_o/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./b006-cys204_cb_ca_c_o-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cys/b006-cys204_cb_ca_c_o/pmf-info.md"

## Visualization

<div id="reduced-view" class="mol-container"></div>
<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('reduced-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/reduced-example.molj", "molj");
    });
});
</script>
