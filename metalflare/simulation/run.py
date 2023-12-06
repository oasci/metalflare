from typing import Any

import argparse
from abc import ABC, abstractclassmethod, abstractmethod
from collections.abc import Collection, Iterable

from loguru import logger

from ..utils import get_obj_from_string
from .contexts import SimulationContextManager


class SimulationRunPrep(ABC):
    r"""Standardized framework for preparing molecular simulation runs."""

    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def prepare_context(
        simulation_context: SimulationContextManager,
    ) -> SimulationContextManager:
        r"""Prepare bash command to run a single simulation. This is unique to each
        simulation package and can be customized based on the context.

        Args:
            simulation_context: Specifies options and parameters for preparing run
                bash command.

        Returns:
            Bash commands in a list to run one stage of a simulation.

        **Examples:**

        TODO:
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_stage_input_lines(cls, context: dict[str, Any]) -> Collection[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        raise NotImplementedError

    @classmethod
    def prepare_slurm_lines(
        cls,
        context: dict[str, Any],
        write: bool = True,
    ) -> list[str]:
        r"""Prepare slurm submission script lines.

        Args:
            context: Specifies options and parameters.

        """
        logger.info("Preparing slurm submission script")
        slurm_lines = ["#!/bin/bash"]
        for k, v in context["sbatch_options"].items():
            slurm_lines.append(f"#SBATCH --{k}={v}")
        slurm_lines.extend(context["slurm_lines"])

        slurm_lines = [l + "\n" for l in slurm_lines if isinstance(l, str)]
        logger.debug("Slurm script:\n{}", "".join(slurm_lines))
        if write:
            logger.info("Writing submissions script at {}", context["slurm_path"])
            with open(context["slurm_path"], mode="w", encoding="utf-8") as f:
                f.writelines(slurm_lines)
        return slurm_lines

    @classmethod
    @abstractmethod
    def get_stage_run_command(cls, context: dict[str, Any]) -> Collection[str]:
        r"""Prepare bash command to run a single simulation.

        Args:
            context: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.

        **Notes:**

        [`prepare_context`][simulation.amber.run.AmberRunPrep.prepare_context]
        should be ran before this.
        """
        raise NotImplementedError

    @classmethod
    @abstractclassmethod
    def prepare_stage(
        cls,
        context: dict[str, Any],
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        raise NotImplementedError

    @classmethod
    @abstractclassmethod
    def prepare(cls, simulation_context: SimulationContextManager) -> None:
        """Run all steps to prepare simulations.

        Args:
            simulation_context: Context manager for simulations.
        """
        raise NotImplementedError


def run_simulation_slurm_prep(
    job_name: str,
    write_dir: str,
    run_path: str,
    slurm_path: str,
    prep_class_string: str,
    simulation_context: SimulationContextManager,
) -> None:
    r"""Run tleap preparation of a system.

    Args:
        job_name: Unique label for these simulations.
        write_dir: Path to local directory where we will write the simulation files.
        run_path: Path to write a bash file that runs all simulations prepared here.
        slurm_path: Path to write slurm submission script.
        prep_class_string: String to a simulation preparation class. For example,
            [`"metalflare.simulation.amber.run.AmberRunPrep"`]
            [simulation.amber.run.AmberRunPrep].
        simulation_context: Context manager for simulations.
    """
    simulation_context.write_dir = write_dir
    simulation_context.slurm_path = slurm_path
    simulation_context.run_path = run_path

    sbatch_options = simulation_context.sbatch_options
    if sbatch_options is not None:
        sbatch_options["job-name"] = job_name
        simulation_context.sbatch_options = sbatch_options
    else:
        raise RuntimeError("sbatch options cannot be None when preparing simulations")

    prep_cls = get_obj_from_string(prep_class_string)  # type: ignore
    prep_cls.prepare_slurm_lines(  # type: ignore
        simulation_context.get(), write=simulation_context.write
    )
    prep_cls.prepare(simulation_context)  # type: ignore


def cli_run_simulation_slurm_prep():
    r"""Command-line interface preparing files for running simulations using slurm."""
    parser = argparse.ArgumentParser(
        description="Prepare files for running simulations using slurm"
    )
    parser.add_argument(
        "job_name",
        type=str,
        nargs="?",
        help="Name of slurm job",
    )
    parser.add_argument(
        "write_dir",
        type=str,
        nargs="?",
        help="Directory to write input files",
    )
    parser.add_argument(
        "run_path",
        type=str,
        nargs="?",
        help="Path to run file",
    )
    parser.add_argument(
        "slurm_path",
        type=str,
        nargs="?",
        help="Path to slurm file",
    )
    parser.add_argument(
        "prep_class_string",
        type=str,
        nargs="?",
        help="String to a simulation preparation class to use.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="*",
        help="Paths to YAML files to use in decreasing precedence.",
    )
    args = parser.parse_args()

    simulation_context = SimulationContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        simulation_context.from_yaml(yaml_path)
    run_simulation_slurm_prep(
        args.job_name,
        args.write_dir,
        args.run_path,
        args.slurm_path,
        args.prep_class_string,
        simulation_context,
    )
