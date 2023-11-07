"""Standardizes residue ID numbering"""
import re

from ..utils import parse_resid, write_resid


def assign_resid(
    line: str, current_resid: int | None, current_original_resid: str
) -> int:
    r"""Determines residue ID based on a consistent numbering scheme.

    Args:
        line: Line that we are determining the residue ID to have.
        current_resid: Current residue ID that we are using.
        current_original_resid: Original residue ID from the PDB file that we are grouping together.

    Returns:
        Assigned residue ID for this line.
    """
    line_resid = parse_resid(line)

    # We have our first residue.
    if current_resid is None:
        try:
            return int(line_resid)
        except ValueError:
            resid = re.sub("[^0-9]", "", str(line_resid))
            if len(str(resid)) == 0:
                resid = "1"
            return int(resid)

    # If the line's residue id is the same as the current original, then we should
    # group this atom with the previous one.
    if line_resid == current_original_resid:
        return current_resid

    # We have the next residue.
    return current_resid + 1


def unify_resid(
    line: str, current_resid: int | None, current_original_resid: str
) -> str:
    r"""Unify residue ID in the PDB line based on previous ones.

    Args:
        line: Line that we are modifying.
        current_resid: Current residue ID that we are using.
        current_original_resid: Original residue ID from the PDB file that we are
            grouping together.

    Returns:
        Line with new residue ID.
    """
    resid = assign_resid(line, current_resid, current_original_resid)
    new_line = write_resid(line, resid)
    return new_line
