# L003: Feature importance

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

## Feature correlation

=== "Reduced"
    <figure markdown>
    ![](./reduced_pls_regression.png){ width=700 }
    </figure>

=== "Oxidized"
    <figure markdown>
    ![](./oxidized_pls_regression.png){ width=700 }
    </figure>

=== "Cu(I)"
    <figure markdown>
    ![](./cu_pls_regression.png){ width=700 }
    </figure>

## Feature importance to ML model

=== "Reduced"
    --8<-- "study/figures/l-his148/l003-his148_hd1-model/reduced-feature-report.md"

=== "Oxidized"
    --8<-- "study/figures/l-his148/l003-his148_hd1-model/oxidized-feature-report.md"

=== "Cu(I)"
    --8<-- "study/figures/l-his148/l003-his148_hd1-model/cu-feature-report.md"
