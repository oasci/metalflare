# b003-cro66_oh-tyr145_hh

TODO:

## Visualization

<div id="b003-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('b003-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/b003.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./b003-cro66_oh-tyr145_hh-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/b-cro-between/b003-cro66_oh-tyr145_hh/pdf-info.md"

### Bandwidth validation

<figure markdown>
![](./b003-cro66_oh-tyr145_hh-hist.svg)
</figure>

## Potential of mean force

<figure markdown>
![](./b003-cro66_oh-tyr145_hh-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/b-cro-between/b003-cro66_oh-tyr145_hh/pmf-info.md"
