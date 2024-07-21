# i002-cro66_oh-tyr145_hh

TODO:

## Visualization

<div id="i002-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('i002-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/i002.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./i002-cro66_oh-tyr145_hh-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/i002-cro66_oh-tyr145_hh/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./i002-cro66_oh-tyr145_hh-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/i002-cro66_oh-tyr145_hh/pmf-info.md"
