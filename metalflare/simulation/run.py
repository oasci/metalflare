from typing import Any

from abc import ABC, abstractmethod
from collections.abc import Iterable

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

    @classmethod
    @abstractmethod
    def get_input_lines(cls, context: dict[str, Any]) -> Iterable[str]:
        r"""Prepare input file lines

        Args:
            context: Specifies options and parameters.
        """

    @classmethod
    @abstractmethod
    def get_stage_run_command(cls, context: dict[str, Any]) -> list[str]:
        r"""Prepare bash command to run a single stage. This is unique to each
        simulation package and can be customized based on the context.

        Args:
            context: Specifies options and parameters.

        Returns:
            Bash commands to run one stage of a simulation.

        **Examples:**

        TODO:
        """