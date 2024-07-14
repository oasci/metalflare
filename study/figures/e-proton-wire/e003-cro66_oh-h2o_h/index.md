# e003-cro66_oh-h2o_h

TODO:

<div id="e003-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('e003-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/e003.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./e003-cro66_oh-h2o_h-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e003-cro66_oh-h2o_h/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e003-cro66_oh-h2o_h-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e003-cro66_oh-h2o_h/pmf-info.md"
