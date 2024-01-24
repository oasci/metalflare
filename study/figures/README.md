# Figures

This folder will contain all figures created during the research process.
It is advised to use vector figures (e.g.,  `.svg` files) when possible, as their size does not change with higher quality outputs.

## Residue indices

Due to our simulation protocol, the numbering system of residues are changed to avoid any mishaps with generating topology files.
In any backend computations or analysis we will use our internal numbering scheme to avoid any inconsistencies.
Frontend communications (e.g., figures) will use the following mapping.

| Backend | Frontend |
| ------- | -------- |
| CRO 65 | CRO 66 |
| TYR 143 | TYR 145 |
| CYM 145 | CYS 147 |
| HID 146 | HIS 148 |
| THR 201 | THR 203 |
| CYM 202 | CYS 204 |
