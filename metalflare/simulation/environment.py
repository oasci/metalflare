from loguru import logger

from .contexts import SimulationContextManager


def get_ion_counts(
    simulation_context: SimulationContextManager,
    system_charge: int | float,
    n_waters: int,
    water_molecule_volume: float = 28.78277638661025,
) -> dict[str, int]:
    r"""Compute the number of cations and anions to achieve desired ionic strength.

    Args:
        simulation_context: A simulation context for system preparation.
        system_charge: Total system charge. Will add counter ions if `neutralize_charge`
            is `True` in `simulation_context`.
        n_waters: Total number of water molecules in the system.
        water_molecule_volume: Approximate volume of a water molecule in cubic
            Angstroms.

    Returns:
        Cations and anions that needs to be added to the system. For example, `2` would
        mean two cations need to be added, but `-2` would mean two anions.
    """
    logger.info("Computing number of extra ions")
    context = simulation_context.get()
    water_box_volume: float = n_waters * water_molecule_volume  # A^3
    logger.debug("Volume of water: {} A^3", water_box_volume)
    water_box_volume /= 1e27  # L
    n_ions = context["solvent_ionic_strength"] * water_box_volume  # moles
    n_ions *= 6.0221409e23  # atoms
    extra_ions = int(round(n_ions, 0))
    ions = {"cations": context["extra_cations"], "anions": context["extra_anions"]}
    ions = {k: v + extra_ions for k, v in ions.items()}
    if context["neutralize_charge"]:
        if system_charge < 0:
            ions["cations"] += abs(int(system_charge))
        elif system_charge > 0:
            ions["anions"] += abs(int(system_charge))
    logger.debug("Ions to add {}", ions)
    return ions
