# G010: Cro66 OH to Thr203 HG1

## Probability density function

<figure markdown>
![](./g010-cro66_oh-thr203_hg1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/g010-cro66_oh-thr203_hg1/pdf-info.md"

<figure markdown>
![](./g010-cro66_oh-thr203_hg1-hist.svg)
</figure>

## Potential of mean force

<figure markdown>
![](./g010-cro66_oh-thr203_hg1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/g010-cro66_oh-thr203_hg1/pmf-info.md"

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
