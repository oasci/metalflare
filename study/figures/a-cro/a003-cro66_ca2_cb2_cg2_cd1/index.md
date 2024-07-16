# a003-cro66_ca2_cb2_cg2_cd1

TODO:

<div id="a003-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('a003-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/a003.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./a003-cro66_ca2_cb2_cg2_cd1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a003-cro66_ca2_cb2_cg2_cd1/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./a003-cro66_ca2_cb2_cg2_cd1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a003-cro66_ca2_cb2_cg2_cd1/pmf-info.md"
