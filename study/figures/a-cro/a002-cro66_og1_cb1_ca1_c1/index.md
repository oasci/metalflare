# A002: Cro66 OG1-CB1-CA1-C1

TODO:

<div id="a001-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('a001-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/a001.molj", "molj");
    });
});
</script>

## Probability density function

<figure markdown>
![](./a002-cro66_og1_cb1_ca1_c1-pdf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a002-cro66_og1_cb1_ca1_c1/pdf-info.md"

## Potential of mean force

TODO:

<figure markdown>
![](./a002-cro66_og1_cb1_ca1_c1-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/a-cro/a002-cro66_og1_cb1_ca1_c1/pmf-info.md"
