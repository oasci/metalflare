# F002: Ser205 -NH to Asn146 =O

TODO:

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
![](./f002-ser205_h-asn146_o-pdf.svg)
</figure>

### Hydrogen bonding

The following table presents the probability of the hydrogen bonding (within 2.5 Å).

| System | H bond |
| ------ | ------ |
| Reduced | 0.047 |
| Oxidized | 0.144 |
| Cu(I) | 0.000 |

### Quantitative

--8<-- "study/figures/f-beta-sheets/f002-ser205_h-asn146_o/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./f002-ser205_h-asn146_o-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/f-beta-sheets/f002-ser205_h-asn146_o/pmf-info.md"