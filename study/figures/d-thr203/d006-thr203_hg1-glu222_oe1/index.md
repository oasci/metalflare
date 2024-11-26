# D006: Thr203 HG1 to Glu222 OE1

## Probability density function

<figure markdown>
![](./d006-thr203_hg1-glu222_oe1-pdf.svg){ width=500 }
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d006-thr203_hg1-glu222_oe1/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./d006-thr203_hg1-glu222_oe1-pmf.svg){ width=500 }
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d006-thr203_hg1-glu222_oe1/pmf-info.md"

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
