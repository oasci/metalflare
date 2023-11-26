import argparse

import MDAnalysis as mda
import numpy as np
from loguru import logger
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from ..structure import get_box_volume


def rotate_positions(positions: np.ndarray, rotation_v: np.ndarray) -> np.ndarray:
    logger.trace("Applying this rotation: {}", rotation_v)
    rotation = Rotation.from_euler("xyz", rotation_v, degrees=True)
    return rotation.apply(positions)


def volume_objective_f(rotation_v: np.ndarray, positions: np.ndarray) -> float:
    r"""Objective function that returns box volume"""
    rotated_positions = rotate_positions(positions, rotation_v)
    return get_box_volume(rotated_positions)


def minimize_box(positions: np.ndarray) -> np.ndarray:
    r""" """
    logger.info("Minimizing box size")
    initial_angles = np.zeros((3,), dtype=np.float64)
    logger.trace("Initial angles: {}", initial_angles)
    results = minimize(
        volume_objective_f,
        initial_angles,
        args=(positions,),
        bounds=[(0, 360)] * 3,
    )
    return results.x


def run_minimize_box(pdb_path: str, output_path: str | None = None) -> np.ndarray:
    r"""Minimize box size by rotating protein.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        Atomic positions rotated to minimum box volume.
    """
    universe = mda.Universe(pdb_path, topology_format="pdb")
    positions = universe.coord.positions
    optimized_rotation_v = minimize_box(positions)
    optimized_positions = rotate_positions(positions, optimized_rotation_v)

    if isinstance(output_path, str):
        universe.coord.positions = optimized_positions
        universe.coord.write(output_path)
    return optimized_positions


def cli_minimize_box() -> None:
    r"""Command-line interface for rotating protein to minimize box volume."""
    parser = argparse.ArgumentParser(
        description="Minimize box size by rotating protein"
    )
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
    run_minimize_box(args.pdb_path, args.output)
