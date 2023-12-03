from typing import Any

import argparse

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from loguru import logger
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from ..structure import get_box_volume


def rotate_positions(
    positions: npt.NDArray[np.float64], rotation_v: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    logger.trace("Applying this rotation: {}", rotation_v)
    rotation: Rotation = Rotation.from_euler("xyz", rotation_v, degrees=True)
    rotated_positions: npt.NDArray[np.float64] = rotation.apply(positions)
    return rotated_positions


def volume_objective_f(
    rotation_v: npt.NDArray[np.float64], positions: npt.NDArray[np.float64]
) -> float:
    r"""Objective function that returns box volume"""
    rotated_positions = rotate_positions(positions, rotation_v)
    return get_box_volume(rotated_positions)


def minimize_box(positions: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    r"""Rotate system to minimize box size.

    Args:
        positions: Cartesian coordinates of system.

    Returns:
        Optimized positions with minimal box size.
    """
    logger.info("Minimizing box size")
    initial_angles = np.zeros((3,), dtype=np.float64)
    logger.trace("Initial angles: {}", initial_angles)
    results = minimize(
        volume_objective_f,
        initial_angles,
        args=(positions,),
        bounds=[(0, 360)] * 3,
    )
    rotated_positions: npt.NDArray[np.float64] = results.x
    return rotated_positions


def run_minimize_box(
    pdb_path: str, output_path: str | None = None
) -> npt.NDArray[np.float64]:
    r"""Minimize box size by rotating protein.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        Atomic positions rotated to minimum box volume.
    """
    u = mda.Universe(pdb_path)
    optimized_rotation_v = minimize_box(u.atoms.positions)
    optimized_positions = rotate_positions(u.atoms.positions, optimized_rotation_v)

    if isinstance(output_path, str):
        u.atoms.positions = optimized_positions
        u.atoms.write(output_path)
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
