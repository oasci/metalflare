"""Standardizes residue ID numbering"""
import argparse
import os
import re
from collections.abc import Iterable

from loguru import logger

from ..utils import parse_resid, replace_in_pdb_line


def assign_resid(
    line: str, current_resid: str | None, prev_original_resid: str
) -> tuple[str, str]:
    r"""Determines residue ID based on a consistent numbering scheme.

    Args:
        line: Line that we are determining the residue ID to have.
        current_resid: Current residue ID that we are using.
        prev_original_resid: Previous residue ID from the PDB file that we are
            grouping together.

    Returns:
        Assigned residue ID for this line.

        Next `prev_original_resid`.
    """
    next_original_resid = parse_resid(line).strip()
    logger.trace("Parsed residue ID from line: {}", next_original_resid)

    # We have our first residue.
    if current_resid is None:
        logger.debug("Current residue ID is None; must be our first atom.")
        assigned_resid = re.sub("[^0-9]", "", str(next_original_resid))
        if len(str(assigned_resid)) == 0:
            assigned_resid = "1"
        logger.trace("Assigning residue ID: {}", assign_resid)
    else:
        # If the line's residue id is the same as the current original, then we should
        # group this atom with the previous one.
        assigned_resid = current_resid
        if next_original_resid != prev_original_resid:
            logger.trace("Parsed residue ID is not the same as previous.")
            logger.trace("Previous residue ID: {}", assigned_resid)
            assigned_resid = str(int(assigned_resid) + 1)
            logger.trace("Next residue ID: {}", assigned_resid)
    return assigned_resid, next_original_resid


def unify_resid(
    line: str, current_resid: str | None, prev_original_resid: str
) -> tuple[str, str, str]:
    r"""Unify residue ID in the PDB line based on previous ones.

    Args:
        line: Line that we are modifying.
        current_resid: Current residue ID that we are using.
        prev_original_resid: Original residue ID from previous line in the PDB file.

    Returns:
        PDB line with the new residue ID.

        Residue ID we are assigning this atom.

        The original residue ID from `line` (i.e., the next `prev_original_resid`).
    """
    assigned_resid, next_original_resid = assign_resid(
        line, current_resid, prev_original_resid
    )
    new_line = replace_in_pdb_line(
        line, next_original_resid, assigned_resid.rjust(4) + " ", 22, 27
    )
    return new_line, assigned_resid, next_original_resid


def run_unify_resids(pdb_path: str, output_path: str | None = None) -> Iterable[str]:
    r"""Unify residue ID in the PDB lines.

    Args:
        pdb_path: Path to PDB file.
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        PDB file lines.
    """
    logger.info("Unify residue IDs from {}", os.path.abspath(pdb_path))
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: list[str] = f.readlines()

    current_resid = None
    current_chain = None
    parse_structure = False
    for i, line in enumerate(pdb_lines):
        logger.trace("Processing line number: {}", i)
        if line.startswith("TER"):
            logger.debug("Encountered 'TER'. Setting parse_structure to False.")
            parse_structure = False
            current_resid = str(parse_resid(pdb_lines[i - 1]))
            continue

        if line.startswith(("ATOM", "HETATM")):
            chain_id = line[21]
            if current_chain is None:
                current_chain = chain_id
            elif current_chain != chain_id:
                current_resid = None
                current_chain = chain_id
            # Activate coordinate parsing on first instance of ATOM or HETATM
            if not parse_structure:
                parse_structure = True
                current_original_resid = str(parse_resid(line))

            pdb_lines[i], current_resid, current_original_resid = unify_resid(
                line, current_resid, current_original_resid
            )

    if output_path is not None:
        logger.info("Writing PDB file to {}", os.path.abspath(output_path))
        with open(output_path, "w+", encoding="utf-8") as f:
            f.writelines(pdb_lines)

    return pdb_lines


def cli_unify_resids() -> None:
    r"""Command-line interface for unifying residue IDs in PDB files."""
    parser = argparse.ArgumentParser(description="Unify residue IDs in PDB")
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
    run_unify_resids(args.pdb_path, args.output)
