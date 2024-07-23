# D004: Thr203 OG1 to H2O

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
![](./d004-thr203_og1-h2o_h-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d004-thr203_og1-h2o_h/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./d004-thr203_og1-h2o_h-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d004-thr203_og1-h2o_h/pmf-info.md"
