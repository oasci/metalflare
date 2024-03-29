# Figures

This folder will contain all figures created during the research process.
It is advised to use vector figures (e.g.,  `.svg` files) when possible, as their size does not change with higher quality outputs.

## Key figures

-   [008-thr203_o-cys204_sg](./008-thr203_o-cys204_sg/)
-   [019-cys204_cb_ca_c_o](./019-cys204_cb_ca_c_o)
-   [012-thr203_hg1_og1_cb_cg2](./012-thr203_hg1_og1_cb_cg2/)
-   [014-ser205_og_cb_ca_n](./014-ser205_og_cb_ca_n/)
-   [011-ser205_hg_og_cb_ca](./011-ser205_hg_og_cb_ca/)
-   [010-cro66_og1_cb1_ca1_c1](./010-cro66_og1_cb1_ca1_c1/)
-   [016-cys204_c-ser205_n_ca_cb](./016-cys204_c-ser205_n_ca_cb/)
-   [017-ser205_h-asn146_o](./017-ser205_h-asn146_o/)

## Residue indices

Due to our simulation protocol, the numbering system of residues are changed to avoid any mishaps with generating topology files.
In any backend computations or analysis we will use our internal numbering scheme to avoid any inconsistencies.
Frontend communications (e.g., figures) will use the following mapping.

| Backend | Frontend |
| ------- | -------- |
| CRO 65 | CRO 66 |
| TYR 143 | TYR 145 |
| ASN 144 | ASN 146 |
| CYM 145 | CYS 147 |
| HID 146 | HIS 148 |
| THR 201 | THR 203 |
| CYM 202 | CYS 204 |
