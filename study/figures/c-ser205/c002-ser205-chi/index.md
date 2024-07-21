# C002: Ser205 χ

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
![](./c002-cys204_c-ser205_n_ca_cb-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/c002-cys204_c-ser205_n_ca_cb/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./c002-cys204_c-ser205_n_ca_cb-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-cys-beta/c002-cys204_c-ser205_n_ca_cb/pmf-info.md"
