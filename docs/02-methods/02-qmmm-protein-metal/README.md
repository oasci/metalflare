# 02-qmmm-protein-metal

## QM/MM molecular dynamics simulation of protein-metal binding in Amber

Add abstract here.

## Prerequisites

While this protocol is generalizable to arbitrary metals and proteins, we designed it with a particular system in mind.
In this section, we define the system-specific properties used as an example.

- `PYTHON_VERSION` = `3.11`
- `CONDA_ENV_NAME` = `metalflare-dev`
- `PDB_ID`: `1JC0`

We often represent these as environment variables that are accessible with `<VAR_NAME>`.

## Virtual environment

A [conda][conda] environment is used to run all code and simulation preparation.
We use a [Makefile][makefile] to automate the creation and management of this environment.
For completeness, we copy the bash commands below.

### Create conda environment

First, we set these environmental variables to simplify the commands.

```bash
export PYTHON_VERSION=3.11 \
export CONDA_NAME=metalflare-dev
```

```bash
conda create -y -n <CONDA_ENV_NAME> python=<PYTHON_VERSION>
```

```bash
conda activate <CONDA_ENV_NAME>
```

### Install dependencies

We use a [conda-lock file][conda-lock file] to specify the exact conda packages installed in the environment.

```bash
conda-lock install -n <CONDA_ENV_NAME> ./conda-lock.yml
```

## Prepare protein structure

### Download PDB file

```bash
wget https://files.rcsb.org/download/<PDB_ID>.pdb
```

[conda]: https://docs.conda.io/en/latest/
[makefile]: https://github.com/oasci/metalflare/blob/main/Makefile
[1JC0]: https://www.rcsb.org/structure/1JC0
[conda-lock file]: https://github.com/oasci/metalflare/blob/main/conda-lock.yml
