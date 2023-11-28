from typing import Any

import os
from collections.abc import Iterable

from loguru import logger

from ..contexts import SimulationContextManager
from ..run import SimulationRunPrep
from .contexts import AmberContextValidator


class AmberRunPrep(SimulationRunPrep):
    r"""Prepares files for Amber simulations"""

    @staticmethod
    def prepare_context(
        simulation_context: SimulationContextManager,
    ) -> SimulationContextManager:
        r"""Preprocessing and validating of context for Amber simulations.

        Args:
            simulation_context: Context manager for simulations.

        """
        logger.info("Preparing context for Amber simulation")
        context = simulation_context.get()

        run_dir = context["output_dir"]

        use_scratch = bool(context["scratch_dir"] is not None)
        if use_scratch:
            run_dir = context["scratch_dir"]
        context["run_dir"] = run_dir

        # Creates all amber-specific paths we will need.
        context["output_path"] = os.path.join(run_dir, context["stage_name"] + ".out")
        context["restart_path"] = os.path.join(run_dir, context["stage_name"] + ".rst")
        context["coord_path"] = os.path.join(run_dir, context["stage_name"] + ".nc")
        context["mdinfo_path"] = os.path.join(
            run_dir, context["stage_name"] + ".mdinfo"
        )

        if use_scratch:
            # We have to use absolute paths with scratch to ensure nothing gets
            # overwritten.
            for k, v in context.items():
                if k[-5:] in ("_path", "_dir"):
                    context[k] = os.path.abspath(v)

        os.makedirs(context["input_dir"], exist_ok=True)
        os.makedirs(context["output_dir"], exist_ok=True)

        simulation_context.update_attributes(context)

        is_valid = AmberContextValidator.validate(simulation_context)
        if not is_valid:
            raise RuntimeError("Context is not valid for Amber!")
        return simulation_context

    @classmethod
    def get_stage_run_command(cls, context: dict[str, Any]) -> list[str]:
        r"""Prepare bash command to run a single Amber simulation using pmemd.

        Args:
            context: Specifies options and parameters.

        Returns:
            Bash commands in a list to run one stage of a simulation.

        **Notes:**

        [`prepare_context`][simulation.amber.run.AmberRunPrep.prepare_context]
        should be ran before this.

        **Uses:**

        The following attributes are possibly used here and should be specified in
        [`simulation_context`][simulation.contexts.SimulationContextManager].

        -   `stage_name`: Unique label for this simulation stage. This will be used
            to build file paths.
        -   `input_path`: Path to input file for this Amber simulations.
        -   `output_dir`: Path to final output directory.
        -   `scratch_dir`: Path to scratch directory if desired. Will run the
            calculations here and then copy to `output_dir`.
        -   `compute_platform`: Computational platform to run the simulation on.
        -   `n_cores`: Number of cores to use for `mpi` simulations if requested.
        """
        use_scratch = bool(context["scratch_dir"] is not None)

        stage_commands = [f"echo 'Starting {context['stage_name']}'", "date"]

        if use_scratch:
            # Adds commands to check if split was already ran.
            check_path = os.path.abspath(
                os.path.join(context["output_dir"], context["stage_name"] + ".rst")
            )
            stage_commands = ["    " + line for line in stage_commands]
            stage_commands.insert(0, f"FILE={check_path}")
            stage_commands.insert(1, 'if [ ! -f "$FILE" ]; then')

        if context["compute_platform"] == "mpi":
            amber_command = f"mpirun -np {context['n_cores']} pmemd.MPI "
        elif context["compute_platform"] == "cuda":
            amber_command = "pmemd.cuda "
        amber_command += f"-O -i {context['input_path']} "
        amber_command += (
            f"-o {context['output_path']} -c {context['prev_restart_path']} "
        )
        amber_command += f"-p {context['topology_path']} "
        amber_command += f"-r {context['restart_path']} -x {context['coord_path']} "
        amber_command += (
            f"-ref {context['ref_coord_path']} -inf {context['mdinfo_path']}"
        )
        logger.debug("Amber command: %s", amber_command[:-2])
        stage_commands.append(amber_command)

        if use_scratch:
            # Add indentation to amber command
            stage_commands[-1] = "    " + stage_commands[-1]

            # Adds commands to move scratch files to output_dir.
            stage_commands.extend(
                [
                    f"    mv {context['output_path']} {context['output_dir']}",
                    f"    mv {context['restart_path']} {context['output_dir']}",
                    f"    mv {context['coord_path']} {context['output_dir']}",
                    "fi",
                ]
            )

        stage_commands.append("")

        return stage_commands

    @classmethod
    def get_stage_input_lines(cls, context: dict[str, Any]) -> list[str]:
        r"""Prepare input file lines for a single stage.

        Args:
            context: Specifies options and parameters.

        Returns:
            Input file lines for a single simulations. The lines do not end in `\n`.
        """
        logger.info("Preparing input lines for stage {}", context["stage_name"])
        input_lines = [context["stage_name"], "&cntrl"]
        for key, value in context["input_kwargs"].items():
            if key in ("restraintmask", "timask1", "scmask1", "timask2", "scmask2"):
                value = f'"{value}"'
            line_to_add = f"    {key}={value},"
            logger.debug("Adding input line: {}", line_to_add.strip())
            input_lines.append(line_to_add)
        input_lines.append("&end")
        return input_lines

    @classmethod
    def prepare_stage(
        cls,
        context: dict[str, Any],
        run_commands: list[str] | None = None,
        write: bool = True,
    ) -> tuple[list[str], list[str]]:
        r"""Write input files for a simulation stage and builds bash commands to
        run all splits.

        Args:
            context: Specifies options and parameters.
            run_commands: Cumulative run commands for all desired stages.
            write: Write input file to disk.

        Returns:
            Input file lines for this stage.

            Updated `run_commands` including this stage.

        **Notes:**

        [`prepare_context`][simulation.amber.run.AmberRunPrep.prepare_context]
        should be ran before this.
        """
        if run_commands is None:
            run_commands = []

        # We do not want to change source context in prepare_context, so we do this
        # here.
        if context["splits"] > 1:
            context["input_kwargs"]["nstlim"] = int(
                context["input_kwargs"]["nstlim"] / context["splits"]
            )

        stage_name = context["stage_name"]
        for i_split in range(1, context["splits"] + 1):
            if context["splits"] > 1:
                stage_name_suffix = f"_split_{i_split:03d}"  # Why n_splits < 1000
            else:
                stage_name_suffix = ""

            context["stage_name"] = stage_name + stage_name_suffix
            stage_input_lines = cls.get_stage_input_lines(context)
            stage_input_path = os.path.join(
                context["input_dir"], context["stage_name"] + ".in"
            )
            if write:
                logger.info("Writing input file at {}", stage_input_path)
                with open(stage_input_path, mode="w", encoding="utf-8") as f:
                    f.writelines([i + "\n" for i in stage_input_lines])

            # Prepare script to run calculations.
            context["coord_path"] = os.path.join(
                context["run_dir"], context["stage_name"] + ".nc"
            )

            logger.info("Adding stage's amber command to run script")
            stage_commands = cls.get_stage_run_command(context)
            run_commands.extend(stage_commands)
        return stage_input_lines, run_commands
