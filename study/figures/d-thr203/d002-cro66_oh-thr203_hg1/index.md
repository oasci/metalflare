# D002: Cro66 OH to Thr203 HG1

TODO:

## Visualization

<div id="d002-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('d002-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/d002.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./d002-cro66_oh-thr203_hg1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/d002-cro66_oh-thr203_hg1/pdf-info.md"

<figure markdown>
![](./d002-cro66_oh-thr203_hg1-hist.svg)
</figure>

## Potential of mean force

<figure markdown>
![](./d002-cro66_oh-thr203_hg1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/d002-cro66_oh-thr203_hg1/pmf-info.md"
