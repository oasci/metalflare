import argparse
from collections.abc import Iterable


def replace_in_pdb_line(
    line: str, orig: str, new: str, start: int | None, stop: int | None
):
    r"""General function to replace parts of a PDB line.

    Args:
        line: PDB line.
        orig: Original value to replace if it exists.
        new: If ``orig`` is in line, replace it with this value.
        start: Slice the line starting here to search.
        stop: Slice the line stopping here to search.
    """
    return line[start:stop].replace(orig, new)


def parse_resid(line: str) -> int:
    r"""Gets the residue ID from a line.

    Args:
        line: Line of a PDB file that starts with ATOM or HETATM.

    Returns:
        Residue ID.
    """
    return int(line[22:30])


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
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_lines: Iterable[str] = f.readlines()

    if record_types is None:
        record_types = ("ATOM", "HETATM", "TER", "END")

    out_lines = keep_lines(pdb_lines, record_types)

    if output_path is not None:
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
