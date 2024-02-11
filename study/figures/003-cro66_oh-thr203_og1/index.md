# 003-cro66_oh-thr203_og1

Looking at threonine (THR) 203 will investigate if the chemical environment around the chromophore changes.

## Visualization

<div id="rogfp-view" class="mol-container"></div>
<script>
var uri = 'https://files.rcsb.org/view/1jc0.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#rogfp-view'),
            { backgroundAlpha: '0.0' }
        );
        let resi1 = 66;
        let atom1Name = "OH";
        let resi2 = 203;
        let atom2Name = "OG1";
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum', opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 66}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 145}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 147}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 148}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 203}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 204}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 205}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'A', resi: 222}, {stick: {}, cartoon: {color: "spectrum", opacity: 0.65}});
        viewer.setStyle({chain: 'B'}, {});
        viewer.setStyle({chain: 'C'}, {});
        viewer.setView([ -180.4790230573923, -10.026709053507197, -44.853784685637685, 111.40397984824831, 0.172578344214806, 0.12220639697970266, 0.9718380201347501, -0.10398641385240419 ]);
        let atom1 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi1, atom: atom1Name}
        )[0];
        let atom2 = viewer.getModel().selectedAtoms(
            {chain: 'A', resi: resi2, atom: atom2Name}
        )[0];
        viewer.addCylinder(
            {
                dashed: true,
                start: {x: atom1.x, y: atom1.y, z: atom1.z},
                end: {x: atom2.x, y: atom2.y, z: atom2.z},
                radius: 0.1,
                color: "#00b4d8"
            }
        );
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

## Potential of mean force

TODO:

<figure markdown>
![](./003-cro66_oh-thr203_og1-pmf.svg)
</figure>

## Probability density function

The figure below shows the distance PDF of CRO66 OH to THR203 OG1.

<figure markdown>
![](./003-cro66_oh-thr203_og1-pdf.svg)
</figure>


rogfp2
y-axis local maxima: [6.64229914e-01 3.36391256e-01 4.07651219e-05]
x-values: [4.03603604 2.73873874 7.85585586]

rogfp2-cu
y-axis local maxima: [1.02555169e+00 5.97934439e-01 9.44273947e-05]
x-values: [2.76576577 3.64864865 5.63963964]

### Bandwidth validation

The figure below shows the histogram and the PDF to validate the bandwidth parameter.
A value of `0.1` was used for both figures.

<figure markdown>
![](./003-cro66_oh-thr203_og1-hist.svg)
</figure>
