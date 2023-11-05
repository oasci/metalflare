# Virtual environment

A [conda][conda] environment is used to run all code and simulation preparation.
We use a [Makefile][makefile] to automate the creation and management of this environment.
For completeness, we copy the bash commands below.

## Create conda environment

```bash
conda create -y -n $(CONDA_ENV_NAME) python=$(PYTHON_VERSION)
```

```bash
conda activate $(CONDA_ENV_NAME)
```

## Install dependencies

We use a [conda-lock file][conda-lock file] to specify the exact conda packages installed in the environment.

```bash
conda-lock install -n $(CONDA_ENV_NAME) ./conda-lock.yml
```

[conda]: https://docs.conda.io/en/latest/
[makefile]: https://github.com/oasci/metalflare/blob/main/Makefile
[conda-lock file]: https://github.com/oasci/metalflare/blob/main/conda-lock.yml
