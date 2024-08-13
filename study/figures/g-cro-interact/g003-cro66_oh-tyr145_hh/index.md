# G003: Cro66 OH - Tyr145 HH

## Probability density function

<figure markdown>
![](./g003-cro66_oh-tyr145_hh-pdf.svg){ width=600 }
</figure>

### Hydrogen bonding

The following table presents the probability of the hydrogen bonding (within 2.5 Å).

| System | H bond |
| ------ | ------ |
| Reduced | 0.612 |
| Oxidized | 0.795 |
| Cu(I) | 0.641 |

### Quantitative

--8<-- "study/figures/g-cro-interact/g003-cro66_oh-tyr145_hh/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./g003-cro66_oh-tyr145_hh-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/g-cro-interact/g003-cro66_oh-tyr145_hh/pmf-info.md"

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
