# Environment

A crucial aspect of consistent development is the standardization of the computational environment.
We use a combination of [`conda`](https://conda.io/) and [`poetry`](https://python-poetry.org/).
Each on its own is more than enough; however, we often want to use packages that are only available in `conda`.
Mixing environment manages like `conda` and `poetry` must be done with care.
This usually involves install all desired `conda` packages and then using only poetry afterwards.
If you want to use a new `conda` package down the road, you normally need to recreate the environment from scratch.

## Steps

### Installing conda

If you do not have `conda` installed, follow the instructions [here](https://docs.conda.io/projects/miniconda/en/latest/#quick-command-line-install).

:::{note}
We recommend using [`libmamba`](https://conda.github.io/conda-libmamba-solver/getting-started/) instead of [`mamba`](https://mamba.readthedocs.io/en/latest/) or [classic `conda`](https://conda.github.io/conda-libmamba-solver/libmamba-vs-classic/).
:::

### Conda environment

First, we setup a `conda` environment called `metalflare-dev`.

```bash
make conda-setup
```

::::{tab-set}

:::{tab-item} From `conda-lock.yml`

```bash
make from-conda-lock
```

:::

:::{tab-item} From scratch

Activate the conda environment.

```bash
conda activate metalflare-dev
```

Add all relevant conda channels so they are exported to `environment.yml`.

```bash
conda config --env --add channels conda-forge
```

Install all desired packages; for example,

```bash
conda install -c conda-forge openmm
```

If needed, write a new `conda-lock` file.

```bash
make write-conda-lock
```

:::

::::

### Poetry-tracked packages

After installing all `conda` packages, we switch over to exclusively using `poetry`.
The following command uses `poetry` to install all packages specified in `pyproject.toml`.

```bash
make install
```

### `pre-commit`

TODO:

```bash
make pre-commit-install
```
