# e008-ser205_og-glh222_he2

TODO:

!!! warning

    Data contained here comes from GLH222 simulations.

!!! quote "Distance"
    <div id="1JC1-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('1JC1-view', {
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
        // viewer.loadStructureFromUrl("/data/005-rogfp-glh-md/structures/protein/1JC0-final.pdb", "pdb");
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/e008.molj", "molj");
    });
});
</script>

## Probability density function

**Bandwidth**: `0.02`.
Anything higher than this starts distorting the relative heights of the peaks.

<figure markdown>
![](./e008-ser205_og-glh222_he2-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e008-ser205_og-glh222_he2/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e008-ser205_og-glh222_he2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e008-ser205_og-glh222_he2/pmf-info.md"
