import MDAnalysis as mda
import numpy as np

from metalflare.pdb.modify.names import run_replace_resnames
from metalflare.pdb.modify.numbering import run_unify_resids
from metalflare.pdb.modify.orientation import minimize_box
from metalflare.pdb.modify.positioning import center_structure
from metalflare.pdb.select import run_select_atoms
from metalflare.pdb.structure import get_box_volume, get_com
from metalflare.pdb.utils import parse_resid, run_filter_pdb


class TestFilterPDBLines:
    def test_1jc0(self, path_1jc0):
        new_pdb = run_filter_pdb(path_1jc0)
        ref_line = "ATOM      1  N   MET A   1     171.453  -1.335  62.904  1.00 58.56           N"
        assert new_pdb[0].strip() == ref_line


class TestRenameResidues:
    def test_1jc0(self, path_1jc0):
        pdb_lines = run_replace_resnames(path_1jc0, {"HOH": "WAT"})
        test_line = pdb_lines[-100]
        assert (
            test_line.strip()
            == "HETATM 5326  O   WAT C 251     167.443  13.783  18.763  1.00 35.94           O"
        )


class TestSelectPDB:
    def test_1jc0_chain(self, path_1jc0):
        atoms = run_select_atoms(path_1jc0, "chainID A")
        assert atoms.positions.shape == (1821, 3)


class TestUnifyPDBResids:
    def test_1jc0(self, path_1jc0):
        # Original is 66
        with open(path_1jc0, "r", encoding="utf-8") as f:
            for line in f:
                if "HETATM" in line and " CRO " in line:
                    resid = int(parse_resid(line.strip()))
                    assert resid == 66

        # Unified is 65
        new_pdb = run_unify_resids(path_1jc0)
        for line in new_pdb:
            if "HETATM" in line and " CRO " in line:
                resid = int(parse_resid(line))
                assert resid == 65


class TestCenterProtein:
    def test_1jc0(self, path_1jc0):
        universe = mda.Universe(path_1jc0, topology_format="pdb")
        atoms = universe.select_atoms("chainID A")
        com = get_com(atoms)
        assert np.allclose(com, np.array([173.95323041, 9.81514191, 43.03364573]))

        centered_atoms = center_structure(atoms)
        com_centered = get_com(centered_atoms)
        assert np.allclose(com_centered, np.array([0.0, 0.0, 0.0]), atol=1e-6)


class TestRotateProtein:
    def test_1jc0(self, path_1jc0):
        universe = mda.Universe(path_1jc0, topology_format="pdb")
        atoms = universe.select_atoms("chainID A")
        positions = atoms.positions
        assert positions.shape == (1821, 3)

        initial_volume = get_box_volume(positions)
        assert initial_volume == 90919.4937923479

        optimized_positions = minimize_box(positions)
        optimized_volume = get_box_volume(optimized_positions)
        assert optimized_volume < initial_volume
