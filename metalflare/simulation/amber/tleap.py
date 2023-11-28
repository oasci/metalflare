r"""Prepare Amber simulations with tleap."""
from typing import Any

import os
import subprocess
import tempfile
from collections.abc import Iterable

from loguru import logger

from ..contexts import SimulationContextManager

FF_WATER_SOLVENT_BOX_MAP: dict[str, Any] = {
    "tip3p": "TIP3PBOX",
    "tip4p": "TIP4PBOX",
    "tip4pew": "TIP4PEWBOX",
    "tip5p": "TIP5PBOX",
    "opc": "OPCBOX",
    "opc3": "OPC3BOX",
    "pol3": "POL3BOX",
    "spce": "SPCBOX",
}
r"""Maps `ff_water` in simulation contexts to tleap box types."""

TLEAP_PATH = os.environ.get("TLEAP_PATH", "tleap")
r"""Path to tleap executable.

You can specify this by setting the path to the `TLEAP_PATH` environmental variable.
For example:

```bash
export TLEAP_PATH="~/miniconda3/envs/metalflare-dev/bin/tleap
```
"""


def get_source_ff_lines(simulation_context: SimulationContextManager) -> list[str]:
    r"""Prepare tleap commands for loading force fields.

    Args:
        simulation_context: A simulation context for system preparation.

    Returns:
        `source leaprc.` commands for tleap.

    **Examples:**

    ```python
    amber_context = {
        "ff_protein": "ff14SB", "ff_water": "tip3p", "ff_small_molecule": "gaff2"
    }
    simulation_context = SimulationContextManager(**amber_context)
    tleap_lines = get_source_ff_lines(simulation_context)
    ```

    would result in

    ```text
    ["source leaprc.protein.ff14SB", "source leaprc.water.tip3p", "source leaprc.gaff2"]
    ```
    """
    logger.info("Creating tleap commands for loading force fields")
    context = simulation_context.get()
    tleap_lines = []
    if isinstance(context["ff_protein"], str):
        tleap_lines.append(f"source leaprc.protein.{context['ff_protein']}")
    if isinstance(context["ff_water"], str):
        tleap_lines.append(f"source leaprc.water.{context['ff_water']}")
    if isinstance(context["ff_dna"], str):
        tleap_lines.append(f"source leaprc.DNA.{context['ff_dna']}")
    if isinstance(context["ff_rna"], str):
        tleap_lines.append(f"source leaprc.RNA.{context['ff_rna']}")
    if isinstance(context["ff_glycam"], str):
        tleap_lines.append(f"source leaprc.{context['ff_glycam']}")
    if isinstance(context["ff_lipid"], str):
        tleap_lines.append(f"source leaprc.{context['ff_lipid']}")
    if isinstance(context["ff_small_molecule"], str):
        tleap_lines.append(f"source leaprc.{context['ff_small_molecule']}")
    if isinstance(context["ff_ions"], str):
        tleap_lines.append(f"source leaprc.{context['ff_ions']}")
    return tleap_lines


def _string_generator(string_list):
    for s in string_list:
        yield s


def parse_tleap_log(log_lines: list[str]) -> dict[str, Any]:
    r"""Parse information from tleap log file.

    Args:
        log_lines: Lines of the `tleap.log` file.

    Returns:
        Parsed information from the log file.

    **Information:**

    -   **`duplicate_atoms`**`: list[dict[str, str]]`

        Each atom in a residue should have a unique type. This list collects information
        of duplicated atom types.

        ```text
        [
            {'residue_number': 23, 'atom_type': 'ND2'},
            {'residue_number': 23, 'atom_type': 'OD1'},
            {'residue_number': 92, 'atom_type': 'NE2'}
        ]
        ```

    """
    logger.info("Parsing tleap log")
    # Initializing information that is not always in tleap.log
    tleap_info: dict[str, Any] = {"duplicate_atoms": [], "unknown_residues": []}
    tleap_generator = _string_generator(log_lines)
    for line in tleap_generator:
        # -- residue 23: duplicate [ ND2] atoms (total 2)
        if "-- residue " in line and "duplicate" in line:
            residue_number = int(line.split()[2][:-1])
            atom_type = line.split("]")[0].split("[")[-1].strip()
            tleap_info["duplicate_atoms"].append(
                {"residue_number": residue_number, "atom_type": atom_type}
            )
        # Unknown residue: CRO   number: 64   type: Nonterminal
        if "Unknown residue:" in line:
            line_split = line.split()
            tleap_info["unknown_residues"].append(
                {
                    "residue_name": line_split[2],
                    "residue_number": line_split[4],
                    "residue_type": line_split[-1],
                }
            )
        # Total vdw box size:                   66.262 69.436 79.756 angstroms.
        # Volume: 366956.734 A^3
        # Mass > 183192.002 amu,  Density > 0.829 g/cc
        # Added 8761 residues.
        if "Total vdw box size: " in line:
            line = next(tleap_generator)
            tleap_info["box_volume"] = float(line.split()[1])
            line = next(tleap_generator)
            tleap_info["box_mass"] = float(line.split()[2])
            tleap_info["box_density"] = float(line.split()[-2])
            while "Added" not in line:
                line = next(tleap_generator)
            tleap_info["n_water_molecules"] = int(line.split()[1])

        # Total unperturbed charge:  -7.000000
        # Total perturbed charge:    -7.000000
        if "Total unperturbed charge:" in line:
            tleap_info["system_charge"] = float(line.strip().split()[-1])

    return tleap_info


def get_prelim_sim_info(
    pdb_path: str,
    simulation_context: SimulationContextManager,
    add_lines: Iterable[str] | None = None,
) -> dict[str, Any]:
    r"""Run a preliminary tleap preparation to get information about the system.

    Args:
        pdb_path: Path to PDB file to load into tleap.
        simulation_context: A simulation context for system preparation.
        add_lines: Additional tleap lines to add before loading the PDB file.

    Returns:
        Parsed information from the tleap log.

    **Examples:**

    The base `tleap` input file is shown below with
    [`AMBER_PROTEIN_STANDARD_CONTEXT`]
    [simulation.amber.contexts.AMBER_PROTEIN_STANDARD_CONTEXT].

    ```bash
    source leaprc.protein.ff19SB
    source leaprc.water.opc3
    <add_lines>
    mol = loadpdb <pdb_path>
    solvatebox mol OPC3BOX 10.0
    savepdb mol <temp file>
    charge mol
    quit
    ```
    """
    context = simulation_context.get()

    # Prepare tleap_lines
    tleap_lines = []
    tleap_lines.extend(get_source_ff_lines(simulation_context))
    if add_lines is None:
        add_lines = []
    tleap_lines.extend(add_lines)
    tleap_lines.append(f"mol = loadpdb {pdb_path}")
    solv_box_name = FF_WATER_SOLVENT_BOX_MAP[context["ff_water"]]
    tleap_lines.append(
        f"solvatebox mol {solv_box_name} {str(context['solvent_padding'])}"
    )
    pdb_output = tempfile.NamedTemporaryFile(mode="r", suffix=".pdb", delete=True)
    tleap_lines.append(f"savepdb mol {pdb_output.name}")
    logger.debug("Setting PDB output to {}", pdb_output.name)
    tleap_lines.extend(["charge mol", "quit"])

    # Writing tleap input file
    tleap_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)
    logger.debug("Writing tleap input to {}", tleap_input.name)
    tleap_string = "\n".join(tleap_lines)
    tleap_input.write(tleap_string)
    logger.debug("tleap input:\n{}", tleap_string)
    tleap_input.close()  # Need to close before running the command.

    logger.info("Running tleap")
    tleap_command = [TLEAP_PATH, "-f", tleap_input.name]
    logger.debug("tleap command: {}", tleap_command)
    completed_process = subprocess.run(
        tleap_command, capture_output=True, text=True, check=False
    )
    os.remove(tleap_input.name)  # Remove temporary input file.
    if completed_process.returncode != 0:
        logger.error("tleap failed!")
    logger.debug("tleap errors:\n{}", completed_process.stderr)
    logger.debug("tleap output:\n{}", completed_process.stdout)

    return parse_tleap_log(completed_process.stdout.split("\n"))


def prepare_amber_files(
    pdb_path: str,
    prmtop_path: str,
    inpcrd_path: str,
    simulation_context: SimulationContextManager,
    pdb_output_path: str | None = None,
    cations: int = 0,
    anions: int = 0,
    add_lines: Iterable[str] | None = None,
) -> dict[str, Any]:
    r"""Prepare amber input files.

    Args:
        pdb_path: Path to PDB file to load into tleap.
        prmtop_path: Path to save topology file.
        inpcrd: Path to save coordinate file.
        simulation_context: A simulation context for system preparation.
        pdb_output_path: Path to save final system.
        cations: Number of cations to add.
        anions: Number of anions to add.
        add_lines: Additional tleap lines to add before generating parameters.

    Returns:
        Parsed information from the tleap log.

    **Examples:**

    The base `tleap` input file is shown below with
    [`AMBER_PROTEIN_STANDARD_CONTEXT`]
    [simulation.amber.contexts.AMBER_PROTEIN_STANDARD_CONTEXT].

    ```bash
    source leaprc.protein.ff19SB
    source leaprc.water.opc3
    <add_lines>
    mol = loadpdb <pdb_path>
    addIons2 mol Na+ <cations>
    addIons2 mol Cl- <anions>
    solvatebox mol OPC3BOX 10.0
    savepdb mol <temp file>
    saveamberparm mol <prmtop_path> <inpcrd_path>
    charge mol
    quit
    ```
    """
    context = simulation_context.get()

    # Prepare tleap_lines
    tleap_lines = []
    tleap_lines.extend(get_source_ff_lines(simulation_context))
    if add_lines is None:
        add_lines = []
    tleap_lines.extend(add_lines)
    tleap_lines.append(f"mol = loadpdb {pdb_path}")
    tleap_lines.extend(
        [
            f"addIons2 mol {context['cation_identity']} {cations}",
            f"addIons2 mol {context['anion_identity']} {anions}",
        ]
    )
    solv_box_name = FF_WATER_SOLVENT_BOX_MAP[context["ff_water"]]
    tleap_lines.append(
        f"solvatebox mol {solv_box_name} {str(context['solvent_padding'])}"
    )
    if pdb_output_path is None:
        pdb_output = tempfile.NamedTemporaryFile(mode="r", suffix=".pdb", delete=True)
        pdb_output_path = pdb_output.name
    tleap_lines.append(f"savepdb mol {pdb_output.name}")
    tleap_lines.append(f"saveamberparm mol {prmtop_path} {inpcrd_path}")
    logger.debug("Setting PDB output to {}", pdb_output.name)
    tleap_lines.extend(["charge mol", "quit"])

    # Writing tleap input file
    tleap_input = tempfile.NamedTemporaryFile(mode="w+", suffix=".in", delete=False)
    logger.debug("Writing tleap input to {}", tleap_input.name)
    tleap_string = "\n".join(tleap_lines)
    tleap_input.write(tleap_string)
    logger.debug("tleap input:\n{}", tleap_string)
    tleap_input.close()  # Need to close before running the command.

    logger.info("Running tleap")
    tleap_command = [TLEAP_PATH, "-f", tleap_input.name]
    logger.debug("tleap command: {}", tleap_command)
    completed_process = subprocess.run(
        tleap_command, capture_output=True, text=True, check=False
    )
    os.remove(tleap_input.name)  # Remove temporary input file.
    if completed_process.returncode != 0:
        logger.error("tleap failed!")
    logger.debug("tleap errors:\n{}", completed_process.stderr)
    logger.debug("tleap output:\n{}", completed_process.stdout)

    return parse_tleap_log(completed_process.stdout.split("\n"))
