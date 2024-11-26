# D007: Thr203 HG1 to Glu222 OE1

## Probability density function

<figure markdown>
![](./d007-thr203_og1-glu222_he2-pdf.svg){ width=500 }
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d007-thr203_og1-glu222_he2/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./d007-thr203_og1-glu222_he2-pmf.svg){ width=500 }
</figure>

### Quantitative

--8<-- "study/figures/d-thr203/d007-thr203_og1-glu222_he2/pmf-info.md"

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
