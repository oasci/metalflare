from metalflare.simulation.environment import get_ion_counts


class TestNumberOfIons:
    def test_1jc0(
        self,
        amber_protein_standard_context,
    ):
        ion_counts = get_ion_counts(
            simulation_context=amber_protein_standard_context,
            system_charge=-7.0,
            n_waters=8761,
        )
        assert ion_counts["cations"] == 30
        assert ion_counts["anions"] == 23
