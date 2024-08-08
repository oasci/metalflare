# J001 Glu222 χ

## Probability density function

<figure markdown>
![](./glu222_n_ca_cb_cg-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/i-tyr145/j001-glu222_n_ca_cb_cg-pdf/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./glu222_n_ca_cb_cg-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/i-tyr145/j001-glu222_n_ca_cb_cg-pdf/pmf-info.md"

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
