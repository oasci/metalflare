# Amber parameters

The following parameters are commonly used in Amber calculations and described in the manual.
We copy them here for convenience.

## `imin`

Flag for running the energy minimization procedure.

-   `= 0` Perform molecular dynamics simulation.
-   `= 1` Perform energy minimization.

## `irest`

Flag to restart a simulation.

-   `= 0` Do not restart the simulation.
-   `= 1` Restart the simulation.

## `ntx`

Option to read the coordinates from `inpcrd` file.

-   `= 1` File is read with no initial velocity information.
-   `= 2` File is read unformatted with no initial velocity information.
-   `= 5` Coordinates and velocities will be read.
    Velocity information will only be used if `irst = 1`.

## `ntmin`

Flag for selecting minimization type.

-   `= 0` Full conjugate gradient minimization.
    The first four cycles are steepest descent at the start of the run and after every nonbonded pair list update.
    Conjugate gradient is slower than steepest descent when far from a minimum. Still, it becomes more efficient when close to the minimum.
-   `= 1` For `ncyc` cycles, the steepest descent method is used then the conjugate gradient is switched on. **[ recommended ]**
-   `= 2` Only the steepest descent method is used.
    This algorithm is fairly popular because it is robust and easy to implement.

## `maxcyc`

Maximum number of minimization cycles allowed.
Values typically range from `1000` to `10000`.
A typical value is `5000`.

## `ncyc`

If `ntmin = 1`, then the minimization method will be switched from steepest descent to conjugate gradient after NCYC cycles.
Values typically range from `50` to `5000`.
A common value is `1000`.

## `ig`

The seed for the pseudo-random number generator.
This affects the starting velocities for MD simulations if `tempi` is nonzero.
If this is `-1`, a random seed will be computed based on the current date and time.
This should almost always be `-1` unless you are reproducing a run.

## `dt`

The time step in picoseconds.
If `ntc = 2` (i.e., SHAKE is on), the maximum value is `0.002` ps (2 fs).
Otherwise, the maximum value is `0.001` ps (1 fs).
Since longer simulations are usually desired, the maximum value is typically used.
Values lower than the maximum can be used, but only if necessary for the phenomena of interest.

## `nstlim`

Number of MD steps to perform.
Multiply `nstlim` and `dt` to get the duration of the MD simulation in picoseconds.

## `nscm`

Flag for the removal of translational and rotational center-of-mass motion every `nscm` steps.
For periodic simulations (`ntb` is `1` or`2`) only the translational center-of-mass motion is removed.
Reasonable values are between `100` and `2000`.
A typical value is `500`.

## `ntr`

Flag for restraining positions of specified atoms using a harmonic potential.

-   `= 0` No constraints.
-   `= 1` Constrain atoms specified in `restraintmask`.

## `restraint_wt`

The weight (in kcal mol <sup>-1</sup> Å<sup>-2</sup>) when `ntr = 1`.
The form of the restraint is $k (\Delta x)^2$ where $\Delta x$ is the deviation of the atom's coordinate from the reference position.
Reasonable values are between `0.1` to `10.0`.
Larger values mean less movement is allowed (i.e., tighter restraint).

## `restraintmask`

Strings that specify the restricted atoms when `ntr = 1`.
To see what atoms will be restrained, you can use `ambmask -p mol.prmtop -c mol.inpcrd -out pdb -find "RESTRAINT_STRING"` in `ambertools`.
Here are some examples of `restraintmask`s and their descriptions.

-   `= "!(@H=)"`: Restrain all atoms except hydrogens.
-   `= "!(:WAT) & !(@H=)"`: Restrain all atoms except hydrogens and water molecules.
-   `= "!(:WAT) & !(@H=) & !(:Na+,Cl-)"`: Same as above, but does not also restrain Na<sup>+</sup> and Cl<sup>-</sup> ions.
-   `= "!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')"`: Restrains protein, DNA, and RNA backbone atoms.
    We must include the `!(:WAT)` to avoid restraining oxygen atoms in the water molecules.
    (TODO: Confirm that atoms after `O` are restraining DNA and RNA backbone atoms.)

## `ntb`

Flab for periodic boundary conditions when computing non-bonded interactions.
Bonds that cross the boundary are not supported.

-   `= 0` No periodic boundary conditions.
-   `= 1` Constant volume.
-   `= 2` Constant pressure.

## `ntf`

Force evaluation type.

-   `= 1` All contributions. **[ required for minimization ]**
-   `= 2` Ignore bond interactions involving hydrogens (when `ntc = 2`).
-   `= 3` All bond interactions are omitted (when `ntc = 3`).

## `ntc`

Flag for SHAKE to perform bond length constraints.

-   `= 1` No SHAKE. **[ required for minimization ]**
-   `= 2` Bonds involving hydrogen are constrained. **[ recommended for MD ]**
-   `= 3` All bonds are constrained.
(Not available for parallel or qmmm runs in sander).

## `cut`

Specifies the nonbonded cutoff in Angstroms.
For particle mesh Ewald (PME), this is used to limit the direct space sum.
`8.0` is usually a good value, but `10.0` is recommended.
(TODO: Why was this value chosen?)
When using an implicit solvent (i.e., `igb > 0`) this should be `9999.0` (effectively infinite).

## `ntt`

Switch for temperature scaling.

-   `= 0` Constant total energy (NVE).
-   `= 1` Constant temperature using the weak-coupling algorithm.
A single scaling factor is used for all atoms.
Generally not recommended.
-   `= 2` Andersen-like temperature coupling scheme, in which imaginary "collisions" are performed with heat bath of temperature `temp0` every `vrand` steps.
-   `= 3` Use Langevin dynamics with the collision frequency `gamma_ln`.
Since Langevin simulations are highly susceptible to "synchronization" artifacts, you should explicitly set `ig` to a different value every restart (e.g., `-1`).
**[ recommended ]**
-   `= 9` Optimized Isokinetic Nose-Hoover chain ensemble (OIN).
Implemented mainly for 3D-RISM and RESPA simulations.
-   `= 10` Stochastic Isokinetic Nose-Hoover RESPA integrator.
Mainly used for RESPA simulations.
-   `= 11` Stochastic version of Berendsen thermostat, also known as the Bussi thermostat.

## `tempi`

Initialization temperature in Kelvin.
If `ntx = 1`, the velocities are assigned from a Maxwellian distribution at `tempi` K.
This has no effect if `ntx = 5`.

## `temp0`

Reference temperature at which the system will be kept in Kelvin.

## `gamma_ln`

The collision frequency, $\gamma$, when `ntt = 3`.
Note that $\gamma$ doesn't need to approximate the physical collision frequency, which is about 50 ps<sup>-1</sup>.
In fact, it is often advantageous to use smaller values around 2 to 5 ps<sup>-1</sup>.
(TODO: More information?)

## `ntp`

Flag for constant pressure dynamics.

-   `= 0` No pressure scaling.
-   `= 1` Isotropic position scaling. **[ recommended ]**
-   `= 2` Anisotropic pressure scaling that can only be used orthogonal boxes.
    Usually, this is only used for membrane simulations, where the surface tensions are different for x, y, and z directions.
    Solutes dissolved in water should not use this.

## `barostat`

Flab used to control which barostat is used in order to control the pressure.

-   `= 1` Berendsen.
-   `= 2` Monte Carlo. **[ recommended ]**

## `pres0`

Reference pressure, in bar, at which the system is maintained.
`1.01325` is almost always used.

## `mcbarint`

Number of steps between volume change attempts performed as part of the Monte Carlo barostat.
`100` is usually the value.

## `comp`

Compressibility of the system when `npt > 0` in units of 10<sup>-6</sup> bar<sup>-1</sup>.
A value of `44.6` is appropriate for water.

## `taup`

Pressure relaxation time in picoseconds when `ntp > 0`.
Recommended values are between `1.0` and `5.0`.
Start with `1.0`, but if your simulations are unstable, you may need to increase this value.

## `ntxo`

Format of the final coordinates, velocities, and box size.

-   `= 1` ASCII.
-   `= 2` Binary NetCDF file. **[ recommended ]**

## `ntwr`

Write the restart file every `ntwr` steps.
A restart file is always written at the end of every run.
One should be cautious about saving too frequently, as the file size will grow quickly, and samples may be time-correlated.
Reasonable values are between `100` and `2000`; however, some situations could call for values outside this range.

## `ntpr`

Print energy information every `ntpr` steps in a human-readable form.
Reasonable values are between `1` and `1000`.

## `ntwx`

Coordinates are written every `ntwx` steps to the `mdcrd` file.
Reasonable values are between `0` and `5000`.
A typical value is `1000`.

## `ioutfm`

Format of coordinate and velocity trajectory files.

-   `= 0` ASCII.
-   `= 1` Binary NetCDF files. **[ recommended ]**

## `iwrap`

Flag for wrapping coordinates around the periodic boundary.

-   `= 0` No wrapping.
-   `= 1` If atoms move across the periodic boundary, amber will wrap them around to the other side. **[ recommended ]**
