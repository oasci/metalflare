from typing import Any

import argparse
import os
from collections.abc import Callable, Iterable

from loguru import logger

from ..utils import parse_atomname, parse_resid, parse_resname, replace_in_pdb_line


def modify_lines(
    pdb_lines: Iterable[str],
    fn_process: Callable[[str, str, str, int | None, int], str],
    fn_args: Iterable[Any],
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    modified_lines = []
    if callable(fn_filter):
        logger.debug("Filter function is provided.")
    for line in pdb_lines:
        if "ATOM" in line or "HETATM" in line:
            if callable(fn_filter):
                filter_return = fn_filter(line).strip()
                logger.trace("Filter function provided: {}", filter_return)
                # Line is in included and should be processed.
                if (include is not None) and (filter_return in include):
                    line = fn_process(line, *fn_args)
                # Line is not in excluded and should be processed.
                if (exclude is not None) and (filter_return not in exclude):
                    line = fn_process(line, *fn_args)
            else:
                line = fn_process(line, *fn_args)
            modified_lines.append(line)

    return modified_lines


def replace_atom_names(
    pdb_lines: Iterable[str], orig_atom_name: str, new_atom_name: str
) -> list[str]:
    r"""Replace all atom names with another.

    Args:
        pdb_lines: List of lines in the PDB file.
        orig_atom_name: Original atom name to replace.
        new_atom_name: New atom name.

    Returns:
        PDB lines with replace atom names.
    """
    orig_atom_name = orig_atom_name.strip().ljust(4)
    new_atom_name = new_atom_name.strip().ljust(4)
    return modify_lines(
        pdb_lines, replace_in_pdb_line, (orig_atom_name, new_atom_name, 13, 17)
    )


def replace_residue_names(
    pdb_lines: Iterable[str],
    orig_resname: str,
    new_resname: str,
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    r"""Replace all instances of residue names with another.

    Args:
        pdb_lines: List of lines in the PDB file.
        orig_resname: Original residue name to replace.
        new_resname: New residue name.
        fn_filter: Filter function to parse part of a line to see if it is in `include`
            or `exclude`.
        include: Only process lines where the result of `fn_filter` is in this list.
        exclude: Do not process lines where the result of `fn_filter` is in this list.

    Returns:
        PDB lines with replace residue names.
    """
    orig_resname = orig_resname.strip().ljust(4)
    new_resname = new_resname.strip().ljust(4)
    logger.info("Renaming '{}' to '{}'", orig_resname, new_resname)
    return modify_lines(
        pdb_lines,
        replace_in_pdb_line,
        (orig_resname, new_resname, 17, 21),
        fn_filter,
        include,
        exclude,
    )


def run_replace_resnames(
    pdb_path: str,
    resname_map: dict[str, str],
    output_path: str | None = None,
    fn_filter: Callable[[str], str] | None = None,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> list[str]:
    r"""Replace residue names.

    Args:
        pdb_path: Path to PDB file.
        resname_map: Original (key) and new (value) mapping of residue names to change.
        output_path: Path to save new PDB file. If `None`, then no file is written.
        fn_filter: Filter function to parse part of a line to see if it is in `include`
            or `exclude`.
        include: Only process lines where the result of `fn_filter` is in this list.
        exclude: Do not process lines where the result of `fn_filter` is in this list.

    Returns:
        PDB lines with changed residues.
    """
    logger.info("Renaming residue names {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    for orig_resname, new_resname in resname_map.items():
        pdb_lines = replace_residue_names(
            pdb_lines, orig_resname, new_resname, fn_filter, include, exclude
        )

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines


def cli_replace_resnames() -> None:
    r"""Command-line interface for renaming residues."""
    parser = argparse.ArgumentParser(description="Rename residues")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "current_resname",
        type=str,
        nargs="?",
        help="Current residue name to replace",
    )
    parser.add_argument(
        "new_resname",
        type=str,
        nargs="?",
        help="New residue name",
    )
    parser.add_argument(
        "--output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.add_argument(
        "--include",
        type=str,
        nargs="*",
        help="Only include these residue indices.",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        nargs="*",
        help="Include all residues except ones with these indices.",
    )
    args = parser.parse_args()
    resname_map = {args.current_resname: args.new_resname}
    if (args.include is not None) or (args.exclude is not None):
        fn_filter = parse_resid
    else:
        fn_filter = None
    run_replace_resnames(
        args.pdb_path, resname_map, args.output, fn_filter, args.include, args.exclude
    )


def run_unify_water_labels(
    pdb_path: str,
    atom_map: dict[str, str] | None = None,
    water_resname: str = "WAT",
    water_atomnames: dict[str, Iterable[str]] | None = None,
    output_path: str | None = None,
) -> Iterable[str]:
    r"""Ensure that water molecule atom names are `O`, `H1`, and `H2`.

    !!! warning

        This has not been tested yet.

    Args:
        pdb_path: Path to PDB file.
        atom_map: Water atom mappings for `O`, `H1`, and `H2`. `H1` and `H2` are the
            first and second hydrogen atoms after `O` in the PDB file, respectively.
            If `None`, then we default to `O`, `H1`, and `H2`.
        water_atomnames: Specifies what water atom names are eligible for replacement.
            If `None`, then we default to `{"O": ["OW"], "H": ["HW"]}`.
        water_resname: Residue name of the water molecules.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        PDB lines with changed residues.
    """
    logger.info("Renaming water atom names in {}", os.path.abspath(pdb_path))

    logger.debug("Water residue name: {}", water_resname)
    if atom_map is None:
        atom_map = {"O": "O", "H1": "H1", "H2": "H2"}
    logger.debug("O atom name: {}", atom_map["O"])
    logger.debug("H1 atom name: {}", atom_map["H1"])
    logger.debug("H2 atom name: {}", atom_map["H2"])
    if water_atomnames is None:
        water_atomnames = {"O": ["OW"], "H": ["HW"]}

    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    for i, line in enumerate(pdb_lines):
        if parse_resname(line).strip() == water_resname:
            original_o_atomname = parse_atomname(line).strip()
            if original_o_atomname in water_atomnames["O"]:
                pdb_lines[i] = replace_atom_names(
                    [line], original_o_atomname, atom_map["O"]
                )[0]

                pdb_lines[i + 1] = replace_atom_names(
                    pdb_lines[i + 1],
                    parse_atomname(pdb_lines[i + 1]).strip(),
                    atom_map["H1"],
                )[0]
                pdb_lines[i + 2] = replace_atom_names(
                    pdb_lines[i + 2],
                    parse_atomname(pdb_lines[i + 2]).strip(),
                    atom_map["H2"],
                )[0]

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines
