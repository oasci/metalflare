import argparse
import os
from collections.abc import Iterable

import MDAnalysis as mda
import numpy as np
from loguru import logger


def replace_in_pdb_line(
    line: str, orig: str, new: str, start: int | None, stop: int | None
) -> str:
    r"""General function to replace parts of a PDB line.

    Args:
        line: PDB line.
        orig: Original value to check if it exists.
        new: If ``orig`` is in line, replace it with this value. This must be formatted
            for all columns, not just the value with no spaces. For example,
            `"   42"` not `"42"`.
        start: Slice the line starting here to replace.
        stop: Slice the line stopping here to replace.
    """
    line_slice = line[start:stop]
    logger.trace("Slice gives us: '{}'", line_slice)
    if orig in line_slice:
        line_slice = new
    return line[:start] + line_slice + line[stop:]


def parse_resid(line: str) -> str:
    r"""Gets the residue ID from a line.

    Args:
        line: Line of a PDB file that starts with ATOM or HETATM.

    Returns:
        Residue ID.
    """
    return line[23:30]


def parse_resname(line: str) -> str:
    r"""Gets the residue name from a line.

    Args:
        line: Line of a PDB file that starts with ATOM or HETATM.

    Returns:
        Residue ID.
    """
    return line[17:21]


def parse_atomname(line: str) -> str:
    r"""Gets the atom name from a line.

    Args:
        line: Line of a PDB file that starts with ATOM or HETATM.

    Returns:
        Residue ID.
    """
    return line[13:17]


def keep_lines(
    lines: Iterable[str],
    record_types: tuple[str, ...] = ("ATOM", "HETATM", "TER", "END"),
) -> list[str]:
    r"""Filter PDB lines to keep in file.

    Args:
        lines: List of lines in the PDB file.
        record_types: Records to keep in the PDB file.

    Returns:
        Filtered lines.
    """
    logger.info("Keeping the following record types: {}", ", ".join(record_types))
    return [line for line in lines if line.startswith(record_types)]


def run_filter_pdb(
    pdb_path: str,
    output_path: str | None = None,
    record_types: tuple[str, ...] | None = None,
) -> list[str]:
    r"""Only keep PDB lines that contain specified record types.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.
        record_types: Records to keep in the PDB file. Defaults to
            `("ATOM", "HETATM", "TER", "END")`.

    Returns:
        PDB file lines.
    """
    logger.info("Filtering PDB lines of {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    if record_types is None:
        record_types = ("ATOM", "HETATM", "TER", "END")
    out_lines = keep_lines(pdb_lines, record_types)

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(out_lines)

    return out_lines


def cli_filter_pdb() -> None:
    r"""Command-line interface for filtering PDB file lines"""
    parser = argparse.ArgumentParser(description="Filter PDB lines")
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
    parser.add_argument(
        "--record-types",
        type=str,
        nargs="*",
        help="Records to keep in the PDB file.",
    )
    args = parser.parse_args()
    run_filter_pdb(args.pdb_path, args.output, args.record_types)


def run_merge_pdbs(*pdb_paths: str, output_path: str | None = None) -> mda.Universe:
    r"""Merge PDB files. No atoms are removed, only added.

    Args:
        *pdb_paths: Paths to PDB files in order of decreasing precedence. We assume
            the residue indices are consistent across PDB files.
    """
    logger.info("Merging PDB files")
    u = mda.Universe(pdb_paths[0])
    atoms_to_add = [mda.Universe(pdb_path).atoms for pdb_path in pdb_paths[1:]]
    u_to_add = mda.core.universe.Merge(*atoms_to_add)

    for residue in u_to_add.residues:
        logger.debug("Processing residue ID of {}", residue.resid)
        available_types = {atom.name for atom in residue.atoms}
        logger.trace("Available atom types: {}", available_types)
        present_types = {atom.name for atom in u.select_atoms(f"resid {residue.resid}")}
        logger.trace("Current atom types: {}", present_types)
        missing_types = {
            atype for atype in available_types if atype not in present_types
        }
        if len(missing_types) > 0:
            logger.trace("Missing atom types: {}", missing_types)
            add_atoms = [atom for atom in residue.atoms if atom.name in missing_types]
            u = mda.core.universe.Merge(u.atoms, mda.AtomGroup(add_atoms))

    # Get the indices of duplicate atoms
    coordinates = u.atoms.positions
    _, unique_indices = np.unique(coordinates, axis=0, return_index=True)
    if len(unique_indices) < len(coordinates):
        logger.info("Cleaning up duplicate atoms")
        duplicate_indices = np.setdiff1d(np.arange(len(coordinates)), unique_indices)
        u.atoms = u.atoms[
            np.isin(np.arange(len(coordinates)), duplicate_indices, invert=True)
        ]

    # Group atoms by residue
    residue_groups = u.atoms.groupby("resids")
    u = mda.core.universe.Merge(
        *[atoms.sort("types") for _, atoms in residue_groups.items()]
    )

    # Remove some Merge artifacts that messes with pdb4amber
    u.del_TopologyAttr("segids")

    if output_path is not None:
        logger.info("Writing merged PDB at {}", output_path)
        u.atoms.write(output_path)
    return u.atoms


def cli_merge_pdbs() -> None:
    r"""Command-line interface for merging PDB files"""
    parser = argparse.ArgumentParser(description="Merge PDB files")
    parser.add_argument(
        "pdb_paths",
        type=str,
        nargs="+",
        help="PDB files to merge",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    args = parser.parse_args()
    if args.output is None:
        raise RuntimeError("--output must be specified")
    run_merge_pdbs(*args.pdb_paths, output_path=args.output)


def run_write_pdb(
    file_paths: Iterable[str],
    output_path: str,
    selection_str: str | None = None,
    stride: int = 1,
) -> None:
    r"""Write PDB file from file paths.

    Args:
        file_paths: Paths of files to load into MDAnalysis.
        output_path: Path to save PDB file.
        selection_str: Selection string for MDAnalysis universe.
    """
    u = mda.Universe(*file_paths)
    atoms = u.select_atoms(selection_str)
    with mda.Writer(output_path, multiframe=True) as W:
        for ts in u.trajectory:
            if ts.frame % stride == 0:
                if isinstance(selection_str, str):
                    atoms = u.select_atoms(selection_str)
                else:
                    atoms = u.atoms
                W.write(atoms)


def cli_write_pdb() -> None:
    r"""Command-line interface for writing PDB from topology and coordinate files."""
    parser = argparse.ArgumentParser(
        description="Write PDB from topology and coordinate files."
    )
    parser.add_argument(
        "output_path",
        type=str,
        nargs="?",
        help="PDB file to write",
    )
    parser.add_argument(
        "--files",
        type=str,
        nargs="+",
        help="Files to load into MDAnalysis universe.",
    )
    parser.add_argument(
        "--selection_str",
        type=str,
        nargs="*",
        help="Selection string for MDAnalysis universe.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        nargs="?",
        help="Stride of trajectory when writing.",
        default=1,
    )
    args = parser.parse_args()
    if args.selection_str is not None:
        args.selection_str = " ".join(args.selection_str)
    run_write_pdb(args.files, args.output_path, args.selection_str, args.stride)
