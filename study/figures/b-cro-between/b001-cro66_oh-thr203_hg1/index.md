# b001-cro66_oh-thr203_hg1

TODO:

## Visualization

<div id="b001-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('b001-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/b001.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./b001-cro66_oh-thr203_hg1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/b001-cro66_oh-thr203_hg1/pdf-info.md"

<figure markdown>
![](./b001-cro66_oh-thr203_hg1-hist.svg)
</figure>

## Potential of mean force

<figure markdown>
![](./b001-cro66_oh-thr203_hg1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/b001-cro66_oh-thr203_hg1/pmf-info.md"
