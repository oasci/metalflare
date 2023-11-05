def coordinate_relevant(line: str) -> bool:
    r"""If a line from a PDB file contains information relevant for coordinates.

    Parameters
    ----------
    line
        String of a PDB file.

    Returns
    -------
        If the line contains atom coordinate information.
    """
    return line.startswith(("ATOM", "HETATM", "END", "TER"))


def parse_resid(line: str) -> int:
    r"""Gets the residue ID from a line.

    Parameters
    ----------
    line
        Line of a PDB file that starts with ATOM or HETATM.

    Returns
    -------
        Residue ID.
    """
    return int(line[22:30])


def set_resid(line: str, resid: int) -> str:
    new_line = line[:22] + str(resid).rjust(4) + " " + line[27:]
    return new_line.strip()
