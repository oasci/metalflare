# i004-tyr145-psi-dist

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

## Probability densities

!!! quote "Reduced"
    <figure markdown>
    ![](./i004-pes-reduced.png)
    </figure>

!!! quote "Oxidized"
    <figure markdown>
    ![](./i004-pes-oxidized.png)
    </figure>

!!! quote "Cu(I)"
    <figure markdown>
    ![](./i004-pes-cu.png)
    </figure>

## Differences

!!! quote "Oxidized vs. reduced"
    <figure markdown>
    ![](./i004-pes-diff-oxd-red.png)
    </figure>

!!! quote "Cu(I) vs. reduced"
    <figure markdown>
    ![](./i004-pes-diff-cu-red.png)
    </figure>
