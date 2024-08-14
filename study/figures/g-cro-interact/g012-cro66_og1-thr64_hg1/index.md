# G012: Cro66 CG2 - Thr64 HG1

## Probability density function

<figure markdown>
![](./g012-cro66_cg2-thr64_hg1-pdf-pdf.svg){ width=600 }
</figure>

### Hydrogen bonding

The following table presents the probability of the hydrogen bonding (within 2.5 Å).

| System | H bond |
| ------ | ------ |
| Reduced | 0.689 |
| Oxidized | 0.977 |
| Cu(I) | 0.895 |
| Na+ | 0.732 |

### Quantitative

--8<-- "study/figures/b-cro-between/b007-cro66_og1-glu222_he2/pdf-info.md"

## Potential of mean force

<figure markdown>
![](./b007-cro66_og1-glu222_he2-pmf.svg)
</figure>

### Quantitative

--8<-- "study/figures/b-cro-between/b007-cro66_og1-glu222_he2/pmf-info.md"

## Visualization

<div id="b003-view" class="mol-container"></div>

<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('b003-view', {
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
        viewer.loadSnapshotFromUrl("/misc/002-molstar-states/b003.molj", "molj");
    });
});
</script>
