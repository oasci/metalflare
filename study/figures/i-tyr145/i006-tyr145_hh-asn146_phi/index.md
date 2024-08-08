# I006: His148-HD1 bonding and Asn146 ψ

## Potential of mean force

=== "Reduced"

    <figure markdown>
    ![](./i006-pes-reduced.png)
    </figure>

=== "Oxidized"

    <figure markdown>
    ![](./i006-pes-oxidized.png)
    </figure>

=== "Cu(I)"

    <figure markdown>
    ![](./i006-pes-cu.png)
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
