#!/usr/bin/env python3

import argparse

from metalflare.simulation.amber.run import AmberRunPrep
from metalflare.simulation.contexts import SimulationContextManager

parser = argparse.ArgumentParser(description="Run tleap")
parser.add_argument(
    "topo_file",
    type=str,
    nargs="?",
    help="Path to .prmtop file",
)
parser.add_argument(
    "coord_file",
    type=str,
    nargs="?",
    help="Path to .inpcrd file",
)
parser.add_argument(
    "slurm_path",
    type=str,
    nargs="?",
    help="Path to slurm file",
)
parser.add_argument(
    "--yaml",
    type=str,
    nargs="*",
    help="Paths to YAML files to use in decreasing precedence.",
)
args = parser.parse_args()

context_manager = SimulationContextManager()
if args.yaml is not None:
    for yaml_path in args.yaml:
        context_manager.from_yaml(yaml_path)

context_manager.slurm_path = args.slurm_path

AmberRunPrep.prepare_slurm_lines(context_manager.get(), write=context_manager.write)
AmberRunPrep.prepare(context_manager)
