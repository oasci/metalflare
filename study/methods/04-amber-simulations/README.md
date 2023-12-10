# 04 - Amber simulations

**Prerequisite(s):** [03-tleap](../03-tleap/)

## Minimization

TODO:

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:1:13"
    ```

### 01_min

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:20:43"
    ```

### 02_min

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:45:67"
    ```

### 03_min

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:69:91"
    ```

### 04_min

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/03-min.yml:93:115"
    ```

## Relation

TODO:

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:1:13"
    ```

### 05_relax_nvt_r

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:19:47"
    ```

### 06_relax_npt_r

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:49:80"
    ```

### 07_relax_npt

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/04-relax.yml:82:111"
    ```

## Production

??? note "YAML"

    ```yaml
    --8<-- "study/methods/04-amber-simulations/05-prod.yml:1:14"
    ```

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
