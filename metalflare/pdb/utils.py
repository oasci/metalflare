import argparse
import os
from collections.abc import Iterable

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
    logger.trace("Slice gives us: `{}`", line_slice)
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


def keep_lines(
    lines: Iterable[str],
    record_types: tuple[str, ...] = ("ATOM", "HETATM", "TER", "END"),
) -> Iterable[str]:
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
    pdb_path: str, output_path: str | None, record_types: tuple[str, ...] | None
) -> Iterable[str]:
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
        pdb_lines: Iterable[str] = f.readlines()

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
