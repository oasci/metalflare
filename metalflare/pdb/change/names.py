from collections.abc import Iterable

from ..utils import replace_in_pdb_line


def modify_lines(pdb_lines, fn_process, fn_args):
    modified_lines = [
        fn_process(line, *fn_args) if "ATOM" in line or "HETATM" in line else line
        for line in pdb_lines
    ]
    return modified_lines


def replace_atom_names(
    pdb_lines: Iterable[str], orig_atom_name: str, new_atom_name: str
) -> Iterable[str]:
    r"""Replace all atom names with another.

    Args:
        pdb_lines: List of lines in the PDB file.
        orig_resname: Original atom name to replace.
        new_resname: New atom name.

    Returns:
        PDB lines with replace atom names.
    """
    orig_atom_name = orig_atom_name.rjust(4)
    new_atom_name = new_atom_name.rjust(4)
    return modify_lines(
        pdb_lines, replace_in_pdb_line, (orig_atom_name, new_atom_name, 13, 17)
    )


def replace_residue_names(
    pdb_lines: Iterable[str], orig_resname: str, new_resname: str
) -> Iterable[str]:
    r"""Replace all instances of residue names with another.

    Args:
        pdb_lines: List of lines in the PDB file.
        orig_resname: Original residue name to replace.
        new_resname: New residue name.

    Returns:
        PDB lines with replace residue names.
    """
    orig_resname = orig_resname.rjust(4)
    new_resname = new_resname.rjust(4)
    return modify_lines(
        pdb_lines, replace_in_pdb_line, (orig_resname, new_resname, 17, 21)
    )
