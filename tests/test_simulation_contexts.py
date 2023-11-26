from metalflare.simulation.amber.contexts import AmberContextValidator


class TestValidateContexts:
    def test_amber_protein_standard_context(self, amber_protein_standard_context):
        AmberContextValidator.validate(amber_protein_standard_context)
