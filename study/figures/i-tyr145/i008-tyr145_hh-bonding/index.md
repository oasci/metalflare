# I008: Tyr145 HH bonding to Cro66-OH vs. Asn146-O

## Potential of mean force

=== "Reduced"

    <figure markdown>
    ![](./pes-reduced.png)
    </figure>

=== "Oxidized"

    <figure markdown>
    ![](./pes-oxidized.png)
    </figure>

=== "Cu(I)"

    <figure markdown>
    ![](./pes-cu.png)
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
