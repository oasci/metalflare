# 03 - tleap

**Prerequisite(s):** [02-protein-prep](../02-protein-prep.md)

[`run_tleap`][simulation.amber.tleap.run_tleap] is used to prepare the topology and coordinate files.

## Force fields

### Protein

TODO:

```yaml
ff_protein: ff19SB
```

### Water

TODO:

```yaml
ff_water: opc3
```

### Ions

TODO:

```yaml
ff_ions: ionslm_126_opc3
```

### Chromophore

The [eGFP chromophore][cro-rcsb] is a non-standard residue that requires external parameters.

<div id="cro-view" class="mol-container"></div>
<script>
var uri = '../../data/001-rogfp-md/structures/protein/1JC0-final.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#cro-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {});
        viewer.setStyle({chain: 'A', resn: 'CRO'}, {stick: {showNonBonded: true}});
        viewer.setView([ -2.0194090909090914, 0.7160909090909091, 0.18300000000000002, 113.69688115360188, 0.4640792686544265, 0.2735166959761386, 0, -0.8425075960651448 ])
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

We use parameters from [this paper][cro-params-paper] that were parameterized with quantum chemical calculations.
They provide [frcmod.xFPchromophores.2022][frcmod.xFPchromophores.2022] and [xFPchromophores.lib.2022][xFPchromophores.lib.2022] that allow us to run classical simulations in Amber.
We use the following commands to load in these files.

```yaml
# Assumes we are in a directory inside of a data directory
add_lines_tleap:
  - 'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
  - '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"} {"ne" "N" "sp2"}'
  - '{"nf" "N" "sp2"} {"ha" "H" "sp3"} {"oh" "O" "sp3"} }'
  - "xFPparams = loadamberparams ../../../methods/03-tleap/cro/frcmod.xFPchromophores.2022"
  - "loadOff ../../../methods/03-tleap/cro/xFPchromophores.lib.2022"

```

## Ions

TODO:

```yaml
cation_identity: Na+
anion_identity: Cl-
neutralize_charge: true
extra_cations: 0
extra_anions: 0
solvent_ionic_strength: 0.150
```

## Solvent

TODO:

```yaml
solvent_padding: 10.0
```

<!-- LINKS -->

[cro-params-paper]: https://doi.org/10.1021/acs.jpcb.3c01486
[frcmod.xFPchromophores.2022]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/methods/03-tleap/cro/frcmod.xFPchromophores.2022
[xFPchromophores.lib.2022]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/methods/03-tleap/cro/xFPchromophores.lib.2022
[cro-rcsb]: https://www.rcsb.org/ligand/CRO
