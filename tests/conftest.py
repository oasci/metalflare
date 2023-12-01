import os

import pytest

from metalflare import enable_logging
from metalflare.simulation.amber.contexts import AMBER_PROTEIN_STANDARD_CONTEXT
from metalflare.simulation.contexts import SimulationContextManager

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture
def test_dir():
    return os.path.abspath(TEST_DIR)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def path_1jc0():
    return os.path.join(TEST_DIR, "files/structures/1JC0.pdb")


@pytest.fixture
def path_1jc0_prepped():
    return os.path.join(TEST_DIR, "files/structures/1JC0-prepped.pdb")


@pytest.fixture
def path_cro_fcrmod():
    return os.path.join(
        TEST_DIR, "files/ff/amber-ff-chromo-params/frcmod.xFPchromophores.2022"
    )


@pytest.fixture
def path_cro_lib():
    return os.path.join(
        TEST_DIR, "files/ff/amber-ff-chromo-params/xFPchromophores.lib.2022"
    )


@pytest.fixture
def amber_protein_standard_context():
    context_manager = SimulationContextManager(**AMBER_PROTEIN_STANDARD_CONTEXT)
    return context_manager


@pytest.fixture
def amber_simulation_standard_context(amber_protein_standard_context):
    context_manager = amber_protein_standard_context
    context_manager.stage_name = "01_min"
    context_manager.input_dir = TEST_DIR
    context_manager.output_dir = TEST_DIR
    context_manager.topo_path = os.path.join(TEST_DIR, "mol.prmtop")
    context_manager.input_path = os.path.join(
        TEST_DIR, f"{context_manager.stage_name}.in"
    )
    context_manager.prev_restart_path = os.path.join(TEST_DIR, "mol.inpcrd")
    context_manager.ref_coord_path = os.path.join(TEST_DIR, "mol.inpcrd")
    context_manager.compute_platform = "mpi"
    context_manager.n_cores = 12
    context_manager.input_kwargs = {
        "irest": 1,
        "ntx": 5,
        "ig": -1,
        "dt": 0.002,
        "nstlim": 500000,
        "nscm": 500,
        "ntr": 1,
        "restraint_wt": 0.5,
        "restraintmask": "!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')",
        "ntb": 2,
        "ntf": 2,
        "ntc": 2,
        "cut": 10.0,
        "ntt": 3,
        "temp0": 300.0,
        "gamma_ln": 5.0,
        "ntp": 1,
        "barostat": 2,
        "pres0": 1.01325,
        "mcbarint": 100,
        "comp": 44.6,
        "taup": 1.0,
        "ntxo": 2,
        "ntwr": 5000,
        "ntpr": 500,
        "ntwx": 5000,
        "ioutfm": 1,
        "iwrap": 1,
    }
    return context_manager
