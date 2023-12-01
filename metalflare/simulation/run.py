from typing import Any

from abc import ABC, abstractclassmethod, abstractmethod
from collections.abc import Collection

from loguru import logger

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
