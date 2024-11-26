import argparse

import MDAnalysis as mda
from loguru import logger

from ..structure import get_com


def center_structure(universe: mda.Universe) -> mda.Universe:
    com = get_com(universe)
    # Shift the entire system so that the COM is at the origin (0, 0, 0)
    logger.info("Centering structure")
    shifted_positions = universe.atoms.positions - com
    universe.atoms.positions = shifted_positions
    return universe


def run_center_structure(pdb_path: str, output_path: str | None = None) -> mda.Universe:
    r"""Center structure by it's center of mass.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        Universe after centering.
    """
    universe = mda.Universe(pdb_path, topology_format="pdb")
    universe = center_structure(universe)

    if isinstance(output_path, str):
        universe.atoms.write(output_path)
    return universe


def cli_center_structure() -> None:
    r"""Command-line interface for centering structure."""
    parser = argparse.ArgumentParser(description="Center structure in PDB file")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    args = parser.parse_args()
    run_center_structure(args.pdb_path, args.output)
