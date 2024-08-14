# E003: H2O O - Cro66 OH

## Water detection

The following table presents the probability of either a water molecule was (a) near or (hydrogen bonding) to the chromophore.

| System | H bond |
| ------ | ------ |
| Reduced | 0.560 |
| Oxidized | 0.601 |
| Cu(I) | 0.538 |
| Na+ | 0.404 |

## Probability density function

<figure markdown>
![](./e003-cro66_oh-h2o_h-pdf.svg){ width=600 }
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e003-cro66_oh-h2o_h/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e003-cro66_oh-h2o_h-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e003-cro66_oh-h2o_h/pmf-info.md"

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
