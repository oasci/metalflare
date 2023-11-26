# Chromophore parameters

These parameters are from [a manuscript](https://pubs.acs.org/doi/10.1021/acs.jpcb.3c01486) that parameterized several chromophores for Amber force fields.

For use with fluorescent protein chromophores.
Tested with EGFP, EBFP, mCherry, ECFP, EYFP, and DsRed proteins
Intended for use with ff14SB, ff19SB and related force fields from the parm94 family.
The chromophores are described using both gaff and parm94 atom types.
The gaff parameters are described by this file and the frcmod and lib files loaded below.
The user must load one of the above parm94-type force fields for those parameters.

Add in gaff atom types that are needed to describe the chromophores.

```bash
addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} {"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }
```

Load in chromophore force field modifications.

```bash
xFPparams = loadamberparams frcmod.xFPchromophores
```

Load in chromophore libraries.

```bash
loadOff xFPchromophores.lib
```
