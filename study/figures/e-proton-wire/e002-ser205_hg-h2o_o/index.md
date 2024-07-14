# e002-ser205_hg-h2o_o

TODO:

<div id="e002-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('e002-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/e002.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./e002-ser205_hg-h2o_o-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e002-ser205_hg-h2o_o/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e002-ser205_hg-h2o_o-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e002-ser205_hg-h2o_o/pmf-info.md"
