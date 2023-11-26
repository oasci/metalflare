import os

import pytest

from metalflare import enable_logging
from metalflare.simulation.amber.contexts import AMBER_PROTEIN_STANDARD_CONTEXT
from metalflare.simulation.contexts import SimulationContextManager

TEST_DIR = os.path.dirname(__file__)


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
    return SimulationContextManager(**AMBER_PROTEIN_STANDARD_CONTEXT)
