# 04 - Amber simulations

**Prerequisite(s):** [03-tleap](../03-tleap/)

The Amber simulation protocol transforms the solvated, parameterized system into a production-ready molecular dynamics trajectory through systematic energy minimization, thermal equilibration, and density relaxation.
This multi-stage approach ensures proper system conditioning while avoiding common simulation artifacts such as structural instabilities, unrealistic conformational changes, or thermodynamic inconsistencies that can compromise the physical validity of the results.

## Simulation Strategy Overview

The complete simulation protocol consists of three distinct phases, each addressing specific aspects of system preparation and equilibration.
The minimization phase removes steric clashes and optimizes unfavorable geometries introduced during system construction.
The relaxation phase gradually brings the system to target temperature and pressure conditions while allowing structural adaptation to the simulation environment.
The production phase collects equilibrium trajectory data under stable thermodynamic conditions for subsequent analysis.
Each phase employs carefully designed restraint strategies, thermodynamic protocols, and sampling parameters to ensure smooth transitions between stages while maintaining structural integrity of critical protein features.
The progressive nature of this approach prevents sudden perturbations that could drive the system into non-physical configurations or cause simulation instabilities.

## Minimization

The energy minimization phase systematically relaxes the initial structure through a series of constrained optimization procedures.
Each minimization stage employs progressively reduced restraints, allowing the system to adapt gradually while preserving essential structural features.
This approach prevents large-scale structural rearrangements that could compromise the integrity of the roGFP2 fold or disrupt critical interactions in the chromophore environment.

### 01_min

The initial minimization stage applies strong positional restraints to all non-hydrogen atoms, allowing optimization of hydrogen positions and relief of immediate steric clashes without significant structural rearrangement.
This conservative approach protects the experimentally determined protein structure while addressing force field inconsistencies and geometric irregularities.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:20:43"
    ```

The restraint weight of 5.0 kcal/mol/Å² provides sufficient constraint to maintain the overall protein architecture while allowing local relaxation.
The restraint mask `!(@H=)` targets all atoms except hydrogens, ensuring that hydrogen positions can optimize freely to establish proper hydrogen bonding networks.
The 5000-step limit with 1000 steepest descent steps followed by conjugate gradient optimization provides adequate convergence for this constrained system.

### 02_min

The second minimization stage relaxes restraints to exclude water molecules and ions, allowing the solvent to optimize around the restrained protein structure.
This selective approach enables proper solvation shell formation while maintaining protein integrity and preventing unrealistic structural distortions.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:45:67"
    ```

The modified restraint mask `!(:WAT) & !(@H=) & !(:Na+,Cl-)` excludes water molecules, hydrogens, and ions from positional restraints while maintaining constraints on protein heavy atoms.
This approach allows the solvent and electrolyte to equilibrate around the protein while preserving the structural integrity of the folded state.

### 03_min

The third minimization stage transitions to backbone-only restraints, allowing side chain relaxation while maintaining secondary structure integrity.
This approach recognizes that side chain conformations may require significant adjustment to accommodate the simulation force field while preserving the essential beta-barrel architecture of the roGFP2 fold.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:69:91"
    ```

The reduced restraint weight of 2.0 kcal/mol/Å² and backbone-specific mask `(@C,CA,N,O,O5',P,O3',C3',C4',C5')` allow greater structural flexibility while maintaining the protein fold.
This restraint pattern targets the polypeptide backbone and nucleic acid structural elements, providing a foundation for side chain optimization.

### 04_min

The final minimization stage further reduces restraint strength while maintaining backbone constraints, allowing near-complete structural relaxation while preserving fold integrity.
This stage provides the optimized starting structure for subsequent thermal equilibration.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:93:115"
    ```

The minimal restraint weight of 1.0 kcal/mol/Å² provides gentle guidance while allowing substantial conformational flexibility.
This final optimization ensures that the system reaches a stable energy minimum before beginning molecular dynamics simulation.

## Relaxation

The equilibration phase gradually brings the minimized system to target simulation conditions through controlled heating and pressure equilibration.
This process requires careful management of temperature and pressure coupling to avoid system instabilities while ensuring proper thermodynamic equilibration.


### 05_relax_nvt_r

The heating phase employs constant volume (NVT) conditions to gradually raise the system temperature from 100 K to 300 K over 20 ps of simulation time.
This controlled heating prevents thermal shock while allowing kinetic energy to distribute throughout the system.
Backbone restraints maintain structural integrity during the heating process.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:19:47"
    ```

The Langevin thermostat (`ntt: 3`) with a collision frequency of 5.0 ps⁻¹ provides efficient temperature control while introducing appropriate stochastic forces.
The initial temperature of 100 K and target temperature of 300 K create a controlled heating ramp over the 10,000-step simulation.
The backbone restraint weight of 1.0 kcal/mol/Å² maintains structural stability during thermal activation.

### 06_relax_npt_r

The pressure equilibration phase transitions to NPT conditions, allowing the system volume to adjust to the target pressure while maintaining temperature control.
Reduced restraints permit greater structural flexibility while ensuring continued fold stability during density equilibration.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:49:80"
    ```

The extended 1 ns simulation (`nstlim: 500000`) allows sufficient time for volume equilibration under the Monte Carlo barostat (`barostat: 2`).
The reduced restraint weight of 0.5 kcal/mol/Å² permits greater conformational flexibility while maintaining essential structural features.
The compressibility of 44.6 × 10⁻⁶ bar⁻¹ approximates the experimental value for liquid water.

### 07_relax_npt

The final equilibration stage removes all positional restraints, allowing the protein to adopt its preferred conformation under the simulation force field and environmental conditions.
This unrestrained equilibration ensures that the system reaches a proper thermodynamic equilibrium before production data collection.

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:82:111"
    ```

The absence of restraints (`ntr: 0`) allows complete conformational freedom, enabling the protein to explore its natural conformational ensemble under simulation conditions.
This 1 ns equilibration provides adequate time for structural relaxation while monitoring system stability and thermodynamic properties.

## Production

The production phase represents the primary data collection period, employing optimized simulation parameters to maximize sampling efficiency while maintaining thermodynamic stability.
The protocol employs GPU acceleration for enhanced computational performance during the extended production trajectories.

The GPU-accelerated configuration employs NVIDIA A100 hardware to achieve optimal performance for long-timescale simulations.
The single-node configuration minimizes communication overhead while providing substantial computational power for production trajectories.

### 08_prod_npt

The simulation employs the NPT ensemble with a time step of `0.002` ps (`2.0` fs) for `50` million MD steps for a total of `100` ns.
Temperature control is achieved through Langevin dynamics with a target temperature of `300` K and a collision frequency of `5.0` ps<sup>-1</sup>.
Pressure control is implemented using the isotropic position scaling method with a Berendsen barostat, targeting a pressure of `1.01325` bar and a pressure relaxation time of `1.0` ps.
Periodic boundary conditions are used for non-bonded interactions with a `10.0` Å cutoff.
Covalent bonds involving hydrogen are constrained with the SHAKE algorithm.
Coordinates are written every `5000` steps (`10` ps) and energy information every `500` steps (`1` ps).

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/05-prod.yml:21:51"
    ```

## Multiple Replica Strategy

The simulation protocol employs three independent replicas to enhance statistical sampling and enable assessment of simulation convergence.
Each replica begins from the same minimized and equilibrated structure but employs different random number seeds, ensuring trajectory divergence and independent sampling of conformational space.

The replica-based approach provides several critical advantages for investigating Cu(I) binding mechanisms.
First, multiple independent trajectories enable assessment of simulation convergence and identification of reproducible structural features.
Second, the enhanced sampling improves statistical precision for calculated ensemble properties such as binding site geometries and hydrogen bonding patterns.
Third, the independent trajectories can reveal rare conformational events or binding site fluctuations that might be missed in single replica simulations.
