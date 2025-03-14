# G006: Cro66 OH to His148 HD1

## Probability density function

<figure markdown>
![](./g006-cro66_oh-his148_hd1-pdf.svg){ width=600 }
</figure>

### Hydrogen bonding

The following table presents the probability of the hydrogen bonding (within 2.5 Å).

| System | H bond |
| ------ | ------ |
| Reduced | 0.486 |
| Oxidized | 0.691 |
| Cu(I) | 0.339 |
| Na+ | 0.347 |

### Quantitative

--8<-- "study/figures/g-cro-interact/g006-cro66_oh-his148_hd1/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./g006-cro66_oh-his148_hd1-pmf.svg){ width=600 }
</figure>

### Quantitative

--8<-- "study/figures/g-cro-interact/g006-cro66_oh-his148_hd1/pmf-info.md"

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
