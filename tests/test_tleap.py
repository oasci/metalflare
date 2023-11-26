import tempfile

from metalflare.simulation.amber.tleap import get_prelim_sim_info, prepare_amber_files


class TestPrelimInfo:
    def test_1jc0(self, path_1jc0_prepped, amber_protein_standard_context):
        tleap_info = get_prelim_sim_info(
            path_1jc0_prepped, amber_protein_standard_context
        )
        assert len(tleap_info["duplicate_atoms"]) == 13
        assert tleap_info["unknown_residues"][-1]["residue_name"] == "CRO"
        assert tleap_info["n_water_molecules"] == 8761
        assert tleap_info["system_charge"] == -6.0

    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_protein_standard_context,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        extra_tleap_lines = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]
        tleap_info = get_prelim_sim_info(
            path_1jc0_prepped,
            amber_protein_standard_context,
            add_lines=extra_tleap_lines,
        )
        assert len(tleap_info["duplicate_atoms"]) == 13
        assert len(tleap_info["unknown_residues"]) == 0
        assert tleap_info["n_water_molecules"] == 8761
        assert tleap_info["system_charge"] == -7.0


class TestPrep:
    def test_1jc0_with_cro(
        self,
        path_1jc0_prepped,
        amber_protein_standard_context,
        path_cro_fcrmod,
        path_cro_lib,
    ):
        extra_tleap_lines = [
            'addAtomTypes { {"cc" "C" "sp2"} {"cd" "C" "sp2"} {"cf" "C" "sp2"} '
            '{"c" "C" "sp2"} {"nd" "N" "sp2"} {"nc" "N" "sp2"}{"ne" "N" "sp2"}'
            '{"nf" "N" "sp2"}{"ha" "H" "sp3"}{"oh" "O" "sp3"} }',
            f"xFPparams = loadamberparams {path_cro_fcrmod}",
            f"loadOff {path_cro_lib}",
        ]
        prmtop_path = tempfile.NamedTemporaryFile().name
        inpcrd_path = tempfile.NamedTemporaryFile().name
        tleap_info = prepare_amber_files(
            path_1jc0_prepped,
            prmtop_path,
            inpcrd_path,
            amber_protein_standard_context,
            cations=30,
            anions=23,
            add_lines=extra_tleap_lines,
        )
        with open(prmtop_path, "r", encoding="utf-8") as f:
            for line in f:
                if r"%FORMAT(10I8) " in line:
                    line = next(f)
                    assert line.split()[0] == "36341"
                    line = next(f)
                    assert line.split()[-2] == "46"
                    line = next(f)
                    assert line.split()[-3] == "1"
                    break
        with open(inpcrd_path, "r", encoding="utf-8") as f:
            for line in f:
                if "default_name" in line:
                    line = next(f)
                    assert line.split()[0] == "36341"
                    line = next(f)
                    assert line.split()[-2] == "26.6466682"
                    break
