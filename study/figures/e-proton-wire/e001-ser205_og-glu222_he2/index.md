# E001: Ser205 OG - Glu222 HE2

## Probability density function

<figure markdown>
![](./e001-ser205_og-glu222_he2-pdf.svg)
</figure>

### Hydrogen bonding

The following table presents the probability of the hydrogen bonding (within 2.5 Å).

| System | H bond |
| ------ | ------ |
| Reduced | 0.299 |
| Oxidized | 0.000 |
| Cu(I) | 0.000 |

### Quantitative

--8<-- "study/figures/e-proton-wire/e001-ser205_og-glu222_he2/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./e001-ser205_og-glu222_he2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/e-proton-wire/e001-ser205_og-glu222_he2/pmf-info.md"

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
