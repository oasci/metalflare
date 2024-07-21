# L002: Cro66 OH to His148 HD1

TODO:

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

## Probability density function

<figure markdown>
![](./l002-cro66_oh-his148_hd1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/l002-cro66_oh-his148_hd1/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./l002-cro66_oh-his148_hd1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/l002-cro66_oh-his148_hd1/pmf-info.md"
