# k002-phe223_ca_c-val224_n_ca

## Probability density function

<figure markdown>
![](./k002-phe223_ca_c-val224_n_ca-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/i-tyr145/k002-phe223_ca_c-val224_n_ca/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./k002-phe223_ca_c-val224_n_ca-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/i-tyr145/k002-phe223_ca_c-val224_n_ca/pmf-info.md"

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
