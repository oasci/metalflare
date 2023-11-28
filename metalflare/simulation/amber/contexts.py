"""Simulation contexts for Amber"""
from typing import Any

from collections.abc import Iterable

from loguru import logger

from ..contexts import ContextValidator

AMBER_PROTEIN_STANDARD_CONTEXT = {
    "ff_protein": "ff19SB",
    "ff_water": "opc3",
    "cation_identity": "Na+",
    "anion_identity": "Cl-",
    "neutralize_charge": True,
    "extra_cations": 0,
    "extra_anions": 0,
    "solvent_ionic_strength": 0.150,
    "solvent_padding": 10.0,
}
r"""Reasonable values for explicitly solvated protein simulations
in water.

The [ff19SB][ff19sb-paper] protein force field offers improved protein backbone
descriptions by incorporating the CMAP approach.
[OPC3][opc3-paper] was found to provide best in-class accuracy and efficiency.

[ff19sb-paper]: https://doi.org/10.1021/acs.jctc.9b00591
[opc3-paper]: https://doi.org/10.1063/1.4960175
"""


# pylint: disable-next=too-few-public-methods
class AmberContextValidator(ContextValidator):
    r"""Validate Amber contexts."""

    ff_protein: Iterable[str] = (
        "ff19SB",
        "ff14SB",
        "ff99SB",
        "ff15ipq",
        "fb15",
        "ff03ua",
    )
    r"""Options for protein force fields."""

    ff_water: Iterable[str] = (
        "tip4p",
        "tip4pew",
        "tip5p",
        "spce",
        "spceb",
        "opc",
        "opc3",
        "opc3pol",
        "opc3pol",
        "pol3",
        "tip3pfb",
        "tip4pfb",
    )
    r"""Options for water force fields."""
    compute_platform: Iterable[str] = ("mpi", "cuda")
    r"""Options for architecture to run simulations on."""

    @staticmethod
    def splits(value: Any, context: dict[str, Any]) -> bool:
        r"""Validate `splits`"""
        if value is not None:
            if value > 1 and context["scratch_dir"] is None:
                logger.error("scratch_dir must be set if splits > 1")
                return False
            if value >= 1000:
                logger.error("splits cannot be larger than 999.")
                return False
        return True
