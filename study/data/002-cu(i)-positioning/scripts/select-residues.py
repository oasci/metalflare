#!/usr/bin/env python3

import argparse

import MDAnalysis as mda

parser = argparse.ArgumentParser(description="Run tleap")
parser.add_argument(
    "topo_path",
    type=str,
    nargs="?",
    help="Topology file",
)
parser.add_argument(
    "coord_path",
    type=str,
    nargs="?",
    help="Coordinate file",
)
parser.add_argument(
    "output",
    type=str,
    nargs="?",
    help="Path to xyz file to write",
)
parser.add_argument(
    "--resi",
    type=str,
    nargs="*",
    help="Residue indices to select",
)
parser.add_argument(
    "--within",
    type=float,
    nargs="?",
    help="Select all residues with this distance in Angstroms.",
    default=5.0,
)
args = parser.parse_args()

u = mda.Universe(args.topo_path, args.coord_path)

residues_select_str = " or ".join(f"resid {i}" for i in args.resi)
selection_atoms = u.select_atoms(residues_select_str)

around_select_str = f"around {args.within} ( " + residues_select_str + " )"
around_atoms = u.select_atoms(around_select_str)

all_atoms = selection_atoms | around_atoms
all_atoms = all_atoms.select_atoms("protein or resname CRO")
all_residues = all_atoms.residues
all_residues.atoms.write(args.output)
