# 03 - tleap

**Prerequisite(s):** [02-protein-prep](../02-protein-prep.md)

The tLeap program serves as the primary tool for constructing complete molecular dynamics systems from prepared protein structures. This stage transforms the cleaned PDB file into a fully parameterized, solvated system ready for energy minimization and molecular dynamics simulation. The [`run_tleap`][simulation.amber.tleap.run_tleap] function orchestrates this complex process through systematic force field assignment, solvation, and ion placement.

## System Construction Overview

The tLeap preparation process addresses several critical aspects of molecular system construction. First, appropriate force field parameters must be assigned to all atoms in the system, including standard amino acids, non-standard residues like the GFP chromophore, and any metal ions or cofactors. Second, the protein must be solvated in an appropriate water model with realistic boundary conditions. Third, counterions must be added to neutralize the system and establish physiological ionic strength. Finally, the complete system topology and initial coordinates must be generated in formats compatible with Amber simulation engines.

The complexity of this process necessitates careful attention to force field compatibility, parameter validation, and system equilibration requirements. Each component of the system requires specific parameterization approaches, and the interactions between different force field components must be consistent and physically reasonable.

## Force Field Selection and Rationale

### Protein Force Field: ff19SB

The ff19SB force field represents the current generation of Amber protein force fields, incorporating significant improvements in backbone and side chain parameterization compared to earlier versions. This force field was specifically designed to address known issues with protein secondary structure stability and folding thermodynamics that affected previous parameter sets.

```yaml
ff_protein: ff19SB
```

The ff19SB parameter set includes refined backbone dihedral parameters derived from high-level quantum mechanical calculations and experimental validation against protein folding data. For roGFP2 systems, this force field provides appropriate treatment of the beta-barrel fold characteristic of GFP-family proteins while maintaining stability of key structural elements including the chromophore environment and metal coordination sites.

The choice of ff19SB is particularly important for studies involving conformational dynamics and long-timescale simulations, as this force field demonstrates improved balance between different secondary structure preferences. This balance is crucial for capturing the subtle conformational changes that may accompany Cu(I) binding and the associated allosteric effects on chromophore environment.

### Water Model: OPC3

The OPC3 water model provides an optimal balance between computational efficiency and physical accuracy for protein simulations. This three-point water model incorporates improved electrostatic parameters and geometry compared to the widely-used TIP3P model, resulting in better reproduction of water's thermodynamic and transport properties.

```yaml
ff_water: opc3
```

OPC3 water demonstrates superior performance in reproducing experimental water density, diffusion coefficients, and dielectric properties across a range of temperatures. For protein simulations, the improved water model leads to more realistic solvation thermodynamics and better representation of water-mediated interactions that can influence protein stability and dynamics.

The selection of OPC3 is particularly relevant for roGFP2 studies because the chromophore environment includes several critical water-mediated hydrogen bonds that stabilize the fluorescent state. The improved electrostatic treatment in OPC3 provides more accurate modeling of these interactions, which may be sensitive to the presence of bound metal ions.

### Ion Parameters: ionslm_126_opc3

Ion parameterization represents one of the most challenging aspects of biomolecular force field development, particularly for transition metals like copper. The ionslm_126_opc3 parameter set provides ion parameters specifically optimized for use with OPC3 water, ensuring consistent thermodynamic properties and appropriate ion-water interactions.

```yaml
ff_ions: ionslm_126_opc3
```

These parameters employ the 12-6-4 Lennard-Jones potential formulation, which includes an additional attractive term to better reproduce ion-water interaction energies and coordination geometries. This approach is particularly important for divalent and transition metal ions, where accurate coordination chemistry is essential for realistic simulation behavior.

For copper ions specifically, these parameters provide reasonable approximations of Cu(I) and Cu(II) behavior in aqueous solution, though the limitations of classical force fields in representing metal coordination chemistry must be acknowledged. The parameters capture the essential electrostatic and size properties of copper ions while maintaining compatibility with the protein and water force field components.

### Chromophore Parameterization: Quantum-Derived Parameters

The eGFP chromophore represents a significant parameterization challenge due to its unique conjugated structure and extended π-electron system. Standard amino acid force fields cannot adequately describe the chromophore's electronic properties, geometric preferences, and electrostatic interactions.

We employ specialized parameters derived from quantum chemical calculations published in [this comprehensive parameterization study][cro-params-paper]. These parameters were developed through systematic quantum mechanical analysis of chromophore conformations, electrostatic properties, and vibrational characteristics, providing a foundation for accurate classical simulation.

The parameterization process involved high-level density functional theory calculations to determine optimal geometry, partial atomic charges, and force constants for all bonded interactions within the chromophore. Particular attention was paid to the conjugated π-system and its interactions with the surrounding protein environment, ensuring proper representation of the electronic structure that underlies the chromophore's optical properties.

```yaml
# Assumes we are in a directory inside of a data directory
add_lines_tleap:
  - 'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
  - '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"} {"ne" "N" "sp2"}'
  - '{"nf" "N" "sp2"} {"ha" "H" "sp3"} {"oh" "O" "sp3"} }'
  - "xFPparams = loadamberparams ../../../methods/03-tleap/cro/frcmod.xFPchromophores.2022"
  - "loadOff ../../../methods/03-tleap/cro/xFPchromophores.lib.2022"
```

The parameter files [frcmod.xFPchromophores.2022][frcmod.xFPchromophores.2022] and [xFPchromophores.lib.2022][xFPchromophores.lib.2022] contain the complete force field definition for the chromophore, including bond, angle, and dihedral parameters as well as partial atomic charges and van der Waals parameters. The atom type definitions extend the standard Amber atom type set to include the specialized hybridization states and electronic environments present in the chromophore.

## Ion Specification and Electrochemical Environment

### Basic Ion Configuration

The ionic composition of the simulation system must accurately reflect physiological conditions while maintaining electrical neutrality and computational stability. The system employs sodium and chloride ions as the primary electrolyte components, representing the dominant ionic species in biological systems.

```yaml
cation_identity: Na+
anion_identity: Cl-
neutralize_charge: true
extra_cations: 0
extra_anions: 0
solvent_ionic_strength: 0.150
```

The neutralization process automatically calculates the net charge of the protein-chromophore system and adds the appropriate number of counterions to achieve overall neutrality. The additional ionic strength of 0.150 M approximates physiological salt concentrations, providing realistic electrostatic screening and ion atmosphere effects.

### Ion Parameter Validation and Selection

Proper ion parameterization requires careful attention to residue and atom naming conventions, which vary between different parameter sets and simulation programs. The ion parameters are defined in the `atomic_ions.lib` file within the Amber installation, and proper identification requires exact matching of residue and atom names.

To determine the appropriate atom and residue names for ions, inspect the `.pixi/envs/dev/dat/leap/lib/atomic_ions.lib` file for entries matching your target ions. For example, the Cu(I) ion parameters are defined as:

```text
!entry.CU1.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "CU" "Cu+" 0 1 131072 1 29 1.000000
```

This specification indicates that Cu(I) ions should be represented with residue name `CU1` and atom name `CU`, with the formal charge of +1. The corresponding PDB entries must match these naming conventions exactly:

```text
HETATM 3944  CU  CU1 A 352    ....
```

Case sensitivity is critical in these specifications, as `Cu1` and `CU1` represent different residue types with potentially different parameterizations. Incorrect naming will result in parameter assignment failures or inappropriate force field treatment.

## Solvation Strategy and Boundary Conditions

### Solvent Box Construction

The solvation process creates a periodic water box around the protein system, establishing appropriate boundary conditions for molecular dynamics simulation. The solvent padding distance determines the minimum separation between protein atoms and box boundaries, ensuring adequate solvation without unnecessary computational overhead.

```yaml
solvent_padding: 10.0
```

A 10.0 Å padding distance provides sufficient solvation to prevent artificial protein-protein interactions across periodic boundaries while maintaining computational efficiency. This distance ensures that the protein experiences a realistic aqueous environment without boundary artifacts that could influence conformational dynamics or metal binding behavior.

The solvation process employs a truncated octahedral box geometry, which provides more efficient packing than rectangular boxes while maintaining proper periodicity. This geometry reduces the total number of water molecules required while ensuring uniform solvation around the protein surface.

## System Validation and Quality Control

### Topology Validation

Following tLeap system construction, comprehensive validation ensures proper parameter assignment and system integrity. The validation process checks for missing parameters, incorrect connectivity, unrealistic geometries, and force field inconsistencies that could compromise simulation stability.

```bash
metalflare-validate-context $YAML_PATH metalflare.simulation.amber.contexts.AmberContextValidator
```

The validation routine examines all bonded interactions to ensure appropriate force constants and equilibrium values, verifies that all atom types have been properly assigned, and checks for common parameterization errors such as missing charges or inappropriate hybridization states.

### Coordinate and Topology Generation

The final step generates Amber-compatible topology and coordinate files that contain all information necessary for molecular dynamics simulation:

```bash
metalflare-tleap $PDB_PATH $TOPO_PATH $COORD_PATH --yaml $YAML_PATH --work_dir "$(dirname "$0")"
```

This process creates the `.prmtop` topology file containing all force field parameters, connectivity information, and atomic properties, as well as the `.inpcrd` coordinate file with initial atomic positions. These files serve as input for all subsequent simulation stages.

### System Verification

A final verification step generates a PDB file from the Amber topology and coordinates, enabling visual inspection of the complete system:

```bash
metalflare-pdb $OUTPUT_PDB_PATH --files $TOPO_PATH $COORD_PATH
```

This verification allows identification of common construction errors such as improper ion placement, incorrect solvation, or parameter assignment failures. Visual inspection can reveal issues that automated validation might miss, such as unusual protein conformations or inappropriate inter-molecular contacts.


<!-- LINKS -->

[cro-params-paper]: https://doi.org/10.1021/acs.jpcb.3c01486
[frcmod.xFPchromophores.2022]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/methods/03-tleap/cro/frcmod.xFPchromophores.2022
[xFPchromophores.lib.2022]: https://gitlab.com/oasci/studies/metalflare/-/blob/main/study/methods/03-tleap/cro/xFPchromophores.lib.2022
[cro-rcsb]: https://www.rcsb.org/ligand/CRO
