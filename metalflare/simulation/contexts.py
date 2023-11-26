from typing import Any

import argparse

from loguru import logger
from ruamel.yaml import YAML

from ..utils import get_obj_from_string


class SimulationContextManager:
    r"""Contexts for setting up molecular simulations."""

    def __init__(self, yaml_path: str | None = None, **kwargs) -> None:
        r"""
        Args:
            yaml_path: Path to YAML file to load into the context.
        """
        # Setting default values
        self.ff_protein: str | None = None
        r"""Molecular mechanics force field used to describe polypeptides."""
        self.ff_water: str | None = None
        r"""Molecular mechanics force fields for water."""
        self.ff_dna: str | None = None
        r"""Molecular mechanics force fields for DNA."""
        self.ff_rna: str | None = None
        r"""Molecular mechanics force fields for RNA."""
        self.ff_glycam: str | None = None
        self.ff_lipid: str | None = None
        r"""Molecular mechanics force fields for lipids."""
        self.ff_small_molecule: str | None = None
        r"""Molecular mechanics force fields for small molecules."""
        self.cation_identity: str = "Na+"
        r"""Identity of the cation in the system."""
        self.anion_identity: str = "Cl-"
        r"""Identity of the anion in the system."""
        self.system_net_charge: int = 0
        r"""Net charge of the molecular system."""
        self.neutralize_charge: bool = True
        r"""Flag to determine if system charge should be neutralized by placing
        additional ions."""
        self.extra_cations: int = 0
        r"""Number of extra cations to add to the system."""
        self.extra_anions: int = 0
        r"""Number of extra anions to add to the system."""
        self.solvent_ionic_strength: float = 0.150
        r"""Ionic strength of the solvent in mole/L."""
        self.solvent_padding: float = 10.0
        r"""Padding between solute and box edge to fill with solvent in Angstroms."""
        self.verbosity: int | str | None = None
        r"""Verbosity level for logging."""
        self.yaml_path = yaml_path
        r"""Path of YAML file that was loaded. Defaults to `None`."""
        self.from_yaml(yaml_path)

        self.update_attributes(kwargs)

    def from_yaml(self, yaml_path: str | None) -> None:
        r"""Load context information from a YAML file. This will only update data
        contained in the YAML file.

        Args:
            yaml_path: Path to YAML file to load.
        """
        if yaml_path is not None:
            logger.info("Loading YAML context from {}", yaml_path)
            yaml = YAML(typ="safe")
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.load(f)
            logger.debug("YAML data:\n{}", yaml_data)
            self.update_attributes(yaml_data)
        self.yaml_path = yaml_path

    def update_attributes(self, attr_dict: dict[str, Any]) -> None:
        r"""Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        for key, value in attr_dict.items():
            setattr(self, key, value)

    def get(self) -> dict[str, Any]:
        r"""Retrieve the context.

        Returns:
            A dictionary representing the current context.
        """
        # The following line filters methods and attributes like __dict__.
        context = {
            k: v for k, v in vars(self).items() if not callable(v) and "__" not in k
        }
        logger.debug("Retrieved context: \n{}", context)
        return context

    def __enter__(self) -> dict[str, Any]:
        r"""Enter the context and return the current context as a dictionary."""
        return self.get()

    def __exit__(self, exc_type, exc_value, exc_tb):
        r"""Exit the context.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback information.
        """


# pylint: disable-next=too-few-public-methods
class ContextValidator:
    r"""Base class for validating simulation contexts."""

    @classmethod
    def validate(cls, context_manager: SimulationContextManager) -> bool:
        r"""Validate contexts for simulations.

        Args:
            context_manager: A simulation context manager to validate.

        Returns:
            If the context is valid.
        """
        logger.info("Validating with: {}", cls.__name__)
        is_valid = True
        for key, value in context_manager.get().items():
            try:
                checker = getattr(cls, key)
            except AttributeError:
                logger.debug("Cannot check {}", key)
                continue
            logger.debug("Checking {}", key)
            if value is not None:
                if isinstance(checker, tuple):
                    if value not in checker:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
                if callable(checker):
                    is_value_valid = checker(value)
                    if not is_value_valid:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
            else:
                logger.debug("  Skipping: None")
        return is_valid

    @staticmethod
    def extra_cations(value: Any) -> bool:
        r"""Validate `extra_cations`"""
        if not isinstance(value, int):
            logger.error("extra_cations must be `int` type")
            return False
        if value < 0:
            logger.error("extra_cations cannot be negative")
            return False
        return True

    @staticmethod
    def extra_anions(value: Any) -> bool:
        r"""Validate `extra_anions`"""
        if not isinstance(value, int):
            logger.error("extra_anions must be `int` type")
            return False
        if value < 0:
            logger.error("extra_anions cannot be negative")
            return False
        return True


def run_context_yaml_validator(yaml_path: str, validator_obj_string: str) -> bool:
    r"""Validate YAML context.

    Args:
        yaml_path: Path to YAML file.
        validator_obj_string: String to validator object.

    Returns:
        If the YAML context is valid.
    """
    logger.info("Validating context from {}", yaml_path)
    validator_cls = get_obj_from_string(validator_obj_string)
    context_manager = SimulationContextManager(yaml_path=yaml_path)
    is_valid = validator_cls.validate(context_manager)
    valid_string = "IS"
    if not is_valid:
        valid_string += " NOT"
    logger.info("Context {} valid", valid_string)
    return is_valid


def cli_validate_yaml_context() -> None:
    r"""Command-line interface for validating YAML context files."""
    parser = argparse.ArgumentParser(description="Validate YAML context file")
    parser.add_argument(
        "yaml_path",
        type=str,
        nargs="?",
        help="Path to YAML file",
    )
    parser.add_argument(
        "validator_string",
        type=str,
        nargs="?",
        help="Import string for validator class",
    )
    args = parser.parse_args()
    is_valid = run_context_yaml_validator(args.yaml_path, args.validator_string)
    if not is_valid:
        raise RuntimeError("Context is not valid")
