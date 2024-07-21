# D001: Thr203 backbone

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

## Potential of mean force

### Reduced

!!! quote ""
    <figure markdown>
    ![](./d001-pes-reduced.png)
    </figure>

### Oxidized

!!! quote ""
    <figure markdown>
    ![](./d001-pes-oxidized.png)
    </figure>

### Cu(I)

!!! quote ""
    <figure markdown>
    ![](./d001-pes-cu.png)
    </figure>
