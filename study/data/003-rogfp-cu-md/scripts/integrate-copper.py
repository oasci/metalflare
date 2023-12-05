#!/usr/bin/env python3

import os

import MDAnalysis as mda

from metalflare.pdb.utils import keep_lines

os.chdir(os.path.dirname(os.path.realpath(__file__)))

TOPO_PATH = "../../001-rogfp-md/simulations/02-prep/mol.prmtop"
COORD_PATH = "../../001-rogfp-md/simulations/02-prep/mol.inpcrd"
CU_XYZ_PATH = "../../002-cu(i)-positioning/calculations/02-dock-copper/xtbopt.xyz"
WRITE_PATH = "../structures/protein/1JCO-cu(i).pdb"

u = mda.Universe(TOPO_PATH, COORD_PATH)
protein = u.select_atoms("protein or resname CRO")
protein.write(WRITE_PATH)

with open(WRITE_PATH, "r", encoding="utf-8") as f:
    pdb_lines = f.readlines()
with open(CU_XYZ_PATH, "r", encoding="utf-8") as f:
    xyz_lines = f.readlines()


def gen_pdb_line(
    atom_type,
    atom_id,
    atom_name,
    residue_name,
    chain_id,
    residue_id,
    x_coord,
    y_coord,
    z_coord,
    element_symbol,
):
    # pylint: disable-next=line-too-long, consider-using-f-string
    line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2s}          {:>2s}{:2s}\n".format(
        str(atom_type),
        atom_id,
        str(atom_name),
        "",
        str(residue_name),
        chain_id,
        residue_id,
        "",
        x_coord,
        y_coord,
        z_coord,
        1.00,
        "",
        str(element_symbol),
        "",
    )
    return line


for line in xyz_lines:
    if "Cu " in line:
        _, x_coord, y_coord, z_coord = line.split()


for i_last_atom in range(len(pdb_lines) - 1, 0, -1):
    if ("ATOM" in pdb_lines[i_last_atom]) or ("HETATM" in pdb_lines[i_last_atom]):
        break
split_line = pdb_lines[i_last_atom].split()
copper_line = gen_pdb_line(
    "HETATM",
    int(split_line[1]) + 1,
    "CU1",
    "CU1",
    "X",
    int(split_line[5]) + 1,
    float(x_coord),
    float(y_coord),
    float(z_coord),
    "Cu",
)
pdb_lines.insert(i_last_atom + 1, copper_line)
for i, line in enumerate(pdb_lines):
    if ("ATOM" in line) or ("HETATM" in line):
        line = list(line)
        line[21] = "A"
        pdb_lines[i] = "".join(line)
pdb_lines = keep_lines(pdb_lines)
with open(WRITE_PATH, "w", encoding="utf-8") as f:
    f.writelines(pdb_lines)
