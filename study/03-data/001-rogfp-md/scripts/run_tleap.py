#!/usr/bin/env python3

import argparse

from metalflare.simulation.amber.tleap import get_prelim_sim_info, prepare_amber_files
from metalflare.simulation.contexts import SimulationContextManager
from metalflare.simulation.environment import get_ion_counts

parser = argparse.ArgumentParser(description="Run tleap")
parser.add_argument(
    "pdb",
    type=str,
    nargs="?",
    help="Path to PDB file",
)
parser.add_argument(
    "topo_path",
    type=str,
    nargs="?",
    help="Where to save topology file",
)
parser.add_argument(
    "coord_path",
    type=str,
    nargs="?",
    help="Where to save coordinate file",
)
parser.add_argument(
    "--yaml",
    type=str,
    nargs="*",
    help="Paths to YAML files to use in decreasing precedence.",
)
args = parser.parse_args()
path_cro_fcrmod = "../ff/cro/frcmod.xFPchromophores.2022"
path_cro_lib = "../ff/cro/xFPchromophores.lib.2022"

context_manager = SimulationContextManager()
if args.yaml is not None:
    for yaml_path in args.yaml:
        context_manager.from_yaml(yaml_path)

extra_tleap_lines = [
    'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
    '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
    '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
    f"xFPparams = loadamberparams {path_cro_fcrmod}",
    f"loadOff {path_cro_lib}",
]

tleap_info = get_prelim_sim_info(
    pdb_path=args.pdb, simulation_context=context_manager, add_lines=extra_tleap_lines
)
ion_counts = get_ion_counts(
    simulation_context=context_manager,
    system_charge=tleap_info["system_charge"],
    n_waters=tleap_info["n_water_molecules"],
)
tleap_info_prep = prepare_amber_files(
    args.pdb,
    prmtop_path=args.topo_path,
    inpcrd_path=args.coord_path,
    simulation_context=context_manager,
    add_lines=extra_tleap_lines,
    cations=ion_counts["cations"],
    anions=ion_counts["anions"],
)
