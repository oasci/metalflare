import argparse

import MDAnalysis as mda


def run_select_atoms(
    pdb_path: str, select_string: str, output_path: str | None = None
) -> mda.core.groups.AtomGroup:
    r"""Select atoms from PDB file.

    Args:
        pdb_path: Path to PDB file.
        select_string: [MDAnalysis selection string]
            (https://docs.mdanalysis.org/stable/documentation_pages/selections.html)
        output_path: Path to save new PDB file. If `None`, then no file is written.

    Returns:
        Atomic positions after atom selection.
    """
    universe = mda.Universe(pdb_path, topology_format="pdb")
    atoms = universe.select_atoms(select_string)

    if isinstance(output_path, str):
        atoms.write(output_path)
    return atoms


def cli_select_atoms() -> None:
    r"""Command-line interface for selecting atoms in PDB file."""
    parser = argparse.ArgumentParser(description="Select atoms from PDB file")
    parser.add_argument(
        "pdb_path",
        type=str,
        nargs="?",
        help="Path to PDB file",
    )
    parser.add_argument(
        "output",
        type=str,
        nargs="?",
        help="Path to new PDB file",
    )
    parser.add_argument(
        "--select_str",
        type=str,
        nargs="+",
        help="Selection string in PDB",
    )
    args = parser.parse_args()
    select_string = " ".join(args.select_str)
    run_select_atoms(args.pdb_path, select_string, args.output)
