# 02-methods

This folder will contain information needed to reproduce experiments, information that would eventually end up in the material and method section of a paper.

It is good practice to keep the information up to date and to enter information as soon as possible.
For instance, reagents shall be best documented at the time of purchase.

## Organization

This is not an exhaustive set of directories.
More can be added (and documented) with numbers `02` to `98`.

### `01-protocols`

TODO: Probably scripts?

### `02-qmmm-protein-metal`

QM/MM molecular dynamic simulation of protein-metal binding in Amber.

### `metalflare`

A modularized package that contains all code used in the project.
This will be installed inside the virtual environment and be callable from scripts throughout the repository.
It's name should change based on this project name.

## Contents

```{toctree}
:glob:

01-protocols/*
*
```
