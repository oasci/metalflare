# G011: Glu222 to Cro66 and Ser205 distance

## Probability densities

=== "Reduced"
    <figure markdown>
    ![](./g011-pes-reduced.png)
    </figure>

=== "Oxidized"
    <figure markdown>
    ![](./g011-pes-oxidized.png)
    </figure>

=== "Cu(I)"
    <figure markdown>
    ![](./g011-pes-cu.png)
    </figure>

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
