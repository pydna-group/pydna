"""The provenance record that ``Source._replay_products`` now runs on."""

from unittest import TestCase

from Bio.Restriction import EcoRI, SalI

from pydna.core.dseqrecord import Dseqrecord
from pydna.methods import (
    cre_lox_excision_or_inversion,
    cre_lox_integration,
    gibson_assembly,
    homologous_recombination_excision_or_inversion,
    homologous_recombination_integration,
    pcr_assembly,
    recombinase_excision_or_inversion,
    recombinase_integration,
    restriction_ligation_assembly,
)
from pydna.methods.cre_lox import LOXP_SEQUENCE
from pydna.methods.recombinase import Recombinase
from pydna.opencloning_models import GibsonAssemblySource, SourceInput
from pydna.primer import Primer
from pydna.provenance import Step, replay
from pydna.provenance.adapters.opencloning import (
    method_name_for,
    source_class_for,
    to_source,
    unsupported,
)

HOMOLOGY_A = "AAGTCCGTTCGTTTTACCTGAAGTCCGTTCGTTTTACCTG"
HOMOLOGY_B = "ATTACAGCATGGGAAGAAAGAATTACAGCATGGGAAGAAA"


def _steps(product):
    return product.source._replay_steps()


class TestStep(TestCase):
    def test_replay_reproduces_the_products(self):
        frag_a = Dseqrecord(HOMOLOGY_A + "cccccc" + HOMOLOGY_B)
        frag_b = Dseqrecord(HOMOLOGY_B + "gggggg" + HOMOLOGY_A)
        product = gibson_assembly([frag_a, frag_b], 40)[0]

        (step,) = _steps(product)
        self.assertEqual(step, Step("gibson", [frag_a, frag_b], {"limit": 40}))
        self.assertEqual(
            [p.seq.seguid() for p in replay(step)],
            [p.seq.seguid() for p in product.source._replay_products()],
        )

    def test_describe(self):
        self.assertEqual(
            Step("gibson", [], {"limit": 40}).describe(), "gibson(limit=40)"
        )


class TestRecoveredParameters(TestCase):
    """Parameters a technique was called with are read back off the record."""

    def test_technique_specific_fields(self):
        backbone = Dseqrecord("cccGAATTCaaaGTCGACccc", circular=True)
        insert = Dseqrecord("ggGAATTCaggtGTCGACgg")
        product = restriction_ligation_assembly(
            [backbone, insert], [EcoRI, SalI], circular_only=True
        )[0]

        (step,) = _steps(product)
        self.assertEqual(step.method, "restriction_ligation")
        self.assertEqual(step.params, {"enzymes": [EcoRI, SalI]})

    def test_pcr_replays_over_template_flanked_by_primers(self):
        template = Dseqrecord("AAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGTTTAAACCCGGG")
        fwd = Primer(str(template.seq)[:20])
        rvs = Primer(str(template.seq.reverse_complement())[:20])
        product = pcr_assembly(template, fwd, rvs, add_primer_features=True)[0]

        (step,) = _steps(product)
        self.assertEqual(step.method, "pcr")
        self.assertEqual(
            [str(s.seq) for s in step.inputs],
            [str(fwd.seq), str(template.seq), str(rvs.seq)],
        )
        self.assertEqual(step.params, {"limit": 20, "add_primer_features": True})


class TestSharedProvenanceClasses(TestCase):
    """OpenCloning records the product, not the technique, for these three."""

    def test_homologous_recombination(self):
        homology = "AAGTCCGTTCGTTTTACCTG"
        genome = Dseqrecord(f"aaaaaa{homology}ccccc{homology}aaaaaa")
        insert = Dseqrecord(f"{homology}gggg{homology}")

        integrated = homologous_recombination_integration(genome, [insert], 20)[0]
        excised = homologous_recombination_excision_or_inversion(genome, 20)[0]

        self.assertEqual(
            [s.method for s in _steps(integrated)],
            ["homologous_recombination_integration"],
        )
        self.assertEqual(
            [s.method for s in _steps(excised)],
            ["homologous_recombination_excision_or_inversion"],
        )

    def test_cre_lox(self):
        genome = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
        insert = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)

        integrated = cre_lox_integration(genome, [insert])[0]
        excised = cre_lox_excision_or_inversion(integrated)[0]

        self.assertEqual(
            [s.method for s in _steps(integrated)], ["cre_lox_integration"]
        )
        self.assertEqual(
            [s.method for s in _steps(excised)], ["cre_lox_excision_or_inversion"]
        )

    def test_recombinase_replays_both_readings_of_a_joined_record(self):
        site1, site2 = "ATGCCCTAAaaCT", "CAaaTTTTTTTCCCT"
        recombinase = Recombinase(site1, site2)
        genome = Dseqrecord("ccccccATGCCCTAAAACTaaaaa")
        insert = Dseqrecord("CAAATTTTTTTCCCTbbbbb", circular=True)

        integrated = recombinase_integration(genome, [insert], recombinase)[0]
        self.assertEqual(
            [s.method for s in _steps(integrated)],
            ["recombinase_integration", "recombinase_assembly"],
        )
        self.assertEqual(
            [s.params for s in _steps(integrated)],
            [{"recombinase": recombinase}] * 2,
        )

        excised = recombinase_excision_or_inversion(
            integrated, recombinase.get_reverse_recombinase()
        )[0]
        self.assertEqual([s.method for s in _steps(excised)], ["recombinase_assembly"])


class TestOpenCloningAdapter(TestCase):
    def test_source_class_for(self):
        self.assertIs(source_class_for("gibson"), GibsonAssemblySource)
        self.assertIsNone(source_class_for("not_a_method"))

    def test_method_name_for(self):
        self.assertEqual(method_name_for(GibsonAssemblySource), "gibson")
        with self.assertRaises(NotImplementedError):
            method_name_for(Dseqrecord)

    def test_to_source(self):
        source = to_source(
            Step("gibson", [], {"limit": 40}),
            input=[SourceInput(sequence=Dseqrecord("ATGC"))],
            circular=True,
        )
        self.assertIsInstance(source, GibsonAssemblySource)
        with self.assertRaises(ValueError):
            to_source(Step("not_a_method"))

    def test_every_registered_method_maps_to_a_source(self):
        self.assertEqual(unsupported(), [])
