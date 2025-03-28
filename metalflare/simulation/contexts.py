from typing import Any

import argparse
from collections.abc import Collection

from loguru import logger
from ruamel.yaml import YAML

from ..utils import get_obj_from_string


class SimulationContextManager:
    r"""Contexts for setting up molecular simulations."""

    # pylint: disable-next=too-many-statements
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
        r"""Molecular mechanics force fields for sugars"""
        self.ff_lipid: str | None = None
        r"""Molecular mechanics force fields for lipids."""
        self.ff_small_molecule: str | None = None
        r"""Molecular mechanics force fields for small molecules."""
        self.ff_ions: str | None = None
        r"""Molecular mechanics force fields for ions."""
        self.cation_identity: str = "Na+"
        r"""Identity of the cation in the system."""
        self.anion_identity: str = "Cl-"
        r"""Identity of the anion in the system."""
        self.system_charge: int = 0
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
        self.scratch_dir: str | None = None
        r"""Specify path for scratch directory if desired. If `None`, we do not use
        scratch."""
        self.compute_platform: str = "mpi"
        r"""Desired platform to run simulations on.

        **Options:**

        -   `mpi`: Message passing interface for central processing units (CPUs).
        -   `cuda`: Compute Unified Device Architecture (CUDA) for graphics processing
            units (GPUs).
        """
        self.cpu_cores: int | None = None
        r"""Number of CPU cores to use if requested"""
        self.input_kwargs: dict[str, Any] | None = None
        r"""Simulation keyword arguments for input files."""
        self.stage_name: str | None = None
        r"""Name or label for simulation stage."""
        self.stages: Collection[dict[str, Any]] | None = None
        """Contexts for successive stages. Stage $i > 0$ is assumed to be restarted from
        stage $i - 1$."""
        self.input_dir: str | None = None
        r"""Path to input directory for current stage."""
        self.input_path: str | None = None
        r"""Path to input file for current stage."""
        self.output_dir: str | None = None
        r"""Path to output directory for current stage."""
        self.output_path: str | None = None
        r"""Path to output file for current stage."""
        self.coord_path: str | None = None
        r"""Path to coordinate file for current stage."""
        self.topo_path: str | None = None
        r"""Path to topology file."""
        self.restart_path: str | None = None
        r"""Path to restart file for this stage."""
        self.prev_coordinate_path: str | None = None
        r"""Path to coordinate file of previous stage"""
        self.prev_restart_path: str | None = None
        r"""Path to restart file from previous stage or initial coordinates."""
        self.ref_coord_path: str | None = None
        r"""Path to reference coordinate file. This is often used for enforcing
        restraints."""
        self.splits: int = 1
        r"""Split simulation stage into several chunks."""
        self.sbatch_options: dict[str, Any] | None = None
        r"""[`sbatch` options](https://slurm.schedmd.com/sbatch.html#SECTION_OPTIONS)
        for a [slurm](https://slurm.schedmd.com/) submission script.
        Some common options are:
        [job-name](https://slurm.schedmd.com/sbatch.html#OPT_job-name),
        [nodes](https://slurm.schedmd.com/sbatch.html#OPT_nodes),
        [ntasks-per-node](https://slurm.schedmd.com/sbatch.html#OPT_ntasks-per-node),
        [cpus-per-task](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-task),
        [gpus](https://slurm.schedmd.com/sbatch.html#OPT_gpus),
        [gres](https://slurm.schedmd.com/sbatch.html#OPT_gres),
        [cpus-per-gpu](https://slurm.schedmd.com/sbatch.html#OPT_cpus-per-gpu),
        [chdir](https://slurm.schedmd.com/sbatch.html#OPT_chdir),
        [output](https://slurm.schedmd.com/sbatch.html#OPT_output),
        [error](https://slurm.schedmd.com/sbatch.html#OPT_error),
        [time](https://slurm.schedmd.com/sbatch.html#OPT_time),
        [clusters](https://slurm.schedmd.com/sbatch.html#OPT_clusters),
        [partition](https://slurm.schedmd.com/sbatch.html#OPT_partition),
        [account](https://slurm.schedmd.com/sbatch.html#OPT_account).

        These options are written in the format of `#SBATCH --{key}={value}`.
        """
        self.slurm_lines: list[str] | None = None
        r"""Lines for a slurm submission script."""
        self.run_path: str | None = None
        r"""Path to run file."""
        self.slurm_path: str | None = None
        r"""Path to slurm submission file."""
        self.write: bool = False
        """Write files."""
        self.write_dir: str | None = None
        """Write directory."""
        self.submit: bool = False
        r"""Submit the job."""
        self.work_dir: str | None = None
        r"""Working directory for preparing calculations."""

        self.yaml_path = yaml_path
        r"""Path of YAML file that was loaded. Defaults to `None`."""
        self.from_yaml(yaml_path)

        self.update(kwargs)

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
            self.update(yaml_data)
        self.yaml_path = yaml_path

    def update(self, attr_dict: dict[str, Any]) -> None:
        r"""Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        logger.debug("Updating context:\n{}", attr_dict)
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
        logger.debug("Retrieved context:\n{}", context)
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

    # pylint: disable=unused-argument

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
        context = context_manager.get()
        for key, value in context.items():
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
                    is_value_valid = checker(value, context)
                    if not is_value_valid:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
            else:
                logger.debug("  Skipping: None")
        return is_valid

    @staticmethod
    def write(value: Any, context: dict[str, Any]) -> bool:
        r"""Validate `write`"""
        if value:
            if context["write_dir"] is None:
                logger.error("write_dir must be set if write is True")
                return False
        return True

    @staticmethod
    def extra_cations(value: Any, context: dict[str, Any]) -> bool:
        r"""Validate `extra_cations`"""
        if not isinstance(value, int):
            logger.error("extra_cations must be `int` type")
            return False
        if value < 0:
            logger.error("extra_cations cannot be negative")
            return False
        return True

    @staticmethod
    def extra_anions(value: Any, context: dict[str, Any]) -> bool:
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
    validator_cls = get_obj_from_string(validator_obj_string)  # type: ignore
    context_manager = SimulationContextManager(yaml_path=yaml_path)
    is_valid: bool = validator_cls.validate(context_manager)  # type: ignore
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
