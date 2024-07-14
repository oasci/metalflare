# e001-ser205_og-glu222_he2

TODO:

<div id="e001-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('e001-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/e001.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./e001-ser205_og-glu222_he2-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e001-ser205_og-glu222_he2/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e001-ser205_og-glu222_he2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e001-ser205_og-glu222_he2/pmf-info.md"
