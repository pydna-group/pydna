from pydna.opencloning_models import SequenceLocationStr
from Bio.SeqFeature import SimpleLocation
from pydna.utils import shift_location
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from unittest import TestCase
from pydna.opencloning_models import (
    AssemblySource,
    AssemblyFragment,
    Source,
    SourceInput,
    TextFileSequence,
    PrimerModel,
    CloningStrategy,
    id_mode,
    get_id,
)
from pydna.primer import Primer
from opencloning_linkml.datamodel import (
    AssemblySource as _AssemblySource,
    AssemblyFragment as _AssemblyFragment,
    Source as _Source,
    SourceInput as _SourceInput,
)
from pydantic import BaseModel, ValidationError
from pydna.assembly2 import (
    golden_gate_assembly,
    crispr_integration,
    ligation_assembly,
    pcr_assembly,
    gateway_assembly,
)
from Bio.Seq import reverse_complement
from Bio.Restriction import BsaI, EcoRI, SalI
import textwrap
import json

# Examples that will be used in several tests ==============================================


# We use this to assign ids to objects in a deterministic way for testing
class Counter:
    def __init__(self, start=1):
        self.counter = start

    def get_next_id(self):
        self.counter += 1
        return str(self.counter)


counter = Counter(start=200)

## Generic assembly examples

example_assembly = [
    (1, SimpleLocation(1, 3), SimpleLocation(4, 6)),
    (2, SimpleLocation(7, 9), SimpleLocation(10, 12)),
]
example_fragments = [
    Dseqrecord("AAAAAAAAAAAAAAAAAAAA"),
    Dseqrecord("TTTTTTTTTTTTTTTTTTT"),
    Primer("AATT", name="forward_primer"),
]

## Golden gate assembly example

insert1 = Dseqrecord("GGTCTCAattaAAAAAttaaAGAGACC", name="insert1")
insert2 = Dseqrecord("GGTCTCAttaaCCCCCatatAGAGACC", name="insert2")
insert3 = Dseqrecord("GGTCTCAatatGGGGGccggAGAGACC", name="insert3")

vector = Dseqrecord("TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT", circular=True, name="vector")

golden_gate_product, *_ = golden_gate_assembly(
    [insert1, insert2, insert3, vector], [BsaI], circular_only=True
)
golden_gate_product.name = "product"

## CRISPR integration example

genome = Dseqrecord("aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac", name="genome")
insert = Dseqrecord("aaccggttAAAAAAAAAttcaaagcac", name="insert")
guide = Primer("ttcaatgcaaacagtaatga", name="guide")

crispr_product, *_ = crispr_integration(genome, [insert], [guide], 8)
crispr_product.name = "product"

## Restriction and ligation example

a = Dseqrecord("aaGAATTCccGTCGACaa")
c, d, e = a.cut(EcoRI, SalI)
ligation_product, *_ = ligation_assembly([c, d, e])

a.name = "a"
a.id = counter.get_next_id()
c.name = "c"
c.id = counter.get_next_id()
d.name = "d"
d.id = counter.get_next_id()
e.name = "e"
e.id = counter.get_next_id()
ligation_product.name = "product"
ligation_product.id = counter.get_next_id()

## PCR

primer1 = Primer("ACGTACGT")
primer2 = Primer(reverse_complement("GCGCGCGC"))

pcr_template = Dseqrecord("ccccACGTACGTAAAAAAGCGCGCGCcccc")

pcr_product, *_ = pcr_assembly(pcr_template, primer1, primer2, limit=8)

primer1.name = "primer1"
primer1.id = counter.get_next_id()
primer2.name = "primer2"
primer2.id = counter.get_next_id()
pcr_template.name = "seq"
pcr_template.id = counter.get_next_id()
pcr_product.name = "product"
pcr_product.id = counter.get_next_id()

## Custom cut

custom_cut_template = Dseqrecord("aaGAATTCccGTCGACaa")
custom_cut_product = custom_cut_template.apply_cut(None, ((3, -4), None))
custom_cut_product.name = "custom_cut_product"

## Gateway

attB1 = "ACAACTTTGTACAAAAAAGCAGAAG"
attP1 = "AAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"

seq1 = Dseqrecord("aaa" + attB1 + "ccc")
seq2 = Dseqrecord("aaa" + attP1 + "ccc")

product_gateway_BP, *_ = gateway_assembly([seq1, seq2], "BP")
product_gateway_BP.name = "product_gateway_BP"

# ========================================================================================


class SequenceLocationStrTest(TestCase):
    def test_methods(self):

        for strand in [1, -1]:
            self.assertEqual(
                SequenceLocationStr.from_biopython_location(
                    SimpleLocation(100, 200, strand)
                ),
                SequenceLocationStr.from_start_and_end(
                    start=100, end=200, strand=strand
                ),
            )
            loc = SequenceLocationStr.from_start_and_end(
                start=20, end=12, strand=strand, seq_len=25
            )
            self.assertEqual(
                loc.to_biopython_location(),
                shift_location(SimpleLocation(20, 37, strand), 0, 25),
            )

        # Test shift
        seq = Dseqrecord("TTTGAGTTGTTTACAACGG", circular=True)
        seq.add_feature(0, 4, type_="CDS", strand=1)
        for shift in range(len(seq)):
            seq_shifted = seq.shifted(shift)
            feat = seq_shifted.features[0]
            loc_pydantic = SequenceLocationStr.from_biopython_location(feat.location)
            loc_biopython = loc_pydantic.to_biopython_location()
            self.assertEqual(loc_biopython, feat.location)

    def test_origin_spanning_start_and_end(self):
        loc = SequenceLocationStr.from_start_and_end(
            start=20, end=10, strand=1, seq_len=30
        )
        self.assertEqual(
            loc.to_biopython_location(),
            shift_location(SimpleLocation(20, 40, strand=1), 0, 30),
        )

        with self.assertRaises(TypeError):
            SequenceLocationStr.from_start_and_end(start=20, end=10, strand=1)

        with self.assertRaises(ValueError):
            SequenceLocationStr.from_start_and_end(
                start=20, end=10, strand=1, seq_len=5
            )

    def test_field_validator(self):
        SequenceLocationStr.field_validator("1..3")
        SequenceLocationStr.field_validator(SequenceLocationStr("1..3"))

        with self.assertRaises(ValueError):
            SequenceLocationStr.field_validator(None)
        with self.assertRaises(ValueError) as e:
            SequenceLocationStr.field_validator(1)
        self.assertEqual(
            e.exception.args[0], "Location must be a string or a SequenceLocationStr"
        )
        with self.assertRaises(ValueError) as e:
            SequenceLocationStr.field_validator("aaa")
        self.assertEqual(e.exception.args[0], "Location 'aaa' is not a valid location")

    def test_pydantic_core_schema(self):
        """Test that the SequenceLocationStr can be used as a field in a Pydantic model."""

        class DummyModel(BaseModel):
            location: SequenceLocationStr

        DummyModel(location="1..3")
        DummyModel(location=SequenceLocationStr("1..3"))


class SourceTest(TestCase):
    def test_to_pydantic_model(self):
        """Test that the Source object is converted to the correct Pydantic model,
        and ids are assigned using the id of the object.
        """
        source = Source(
            input=[
                SourceInput(sequence=example_fragments[0]),  # Dseqrecord
                AssemblyFragment(
                    sequence=example_fragments[1],
                    left_location=None,
                    right_location=SimpleLocation(10, 12),
                    reverse_complemented=True,
                ),
                SourceInput(sequence=example_fragments[2]),  # Primer
            ],
        )
        self.assertEqual(
            source.to_pydantic_model(1),
            _Source(
                id=1,
                input=[
                    _SourceInput(sequence=id(example_fragments[0])),
                    _AssemblyFragment(
                        sequence=id(example_fragments[1]),
                        left_location=None,
                        right_location="11..12",
                        reverse_complemented=True,
                    ),
                    _SourceInput(sequence=id(example_fragments[2])),
                ],
            ),
        )

    def test_history_string(self):
        # If source is None, returns empty string
        self.assertEqual(Dseqrecord("A").history(), "")
        self.assertEqual(
            golden_gate_product.history(),
            textwrap.dedent(
                """
            ╙── product (Dseqrecord(o39))
                └─╼ RestrictionAndLigationSource
                    ├─╼ insert1 (Dseqrecord(-27))
                    ├─╼ insert2 (Dseqrecord(-27))
                    ├─╼ insert3 (Dseqrecord(-27))
                    └─╼ vector (Dseqrecord(o35))
            """
            ).strip(),
        )

        # Handles primers as well
        self.assertEqual(
            crispr_product.history(),
            textwrap.dedent(
                """
            ╙── product (Dseqrecord(-27))
                └─╼ CRISPRSource
                    ├─╼ genome (Dseqrecord(-45))
                    ├─╼ insert (Dseqrecord(-27))
                    └─╼ guide (id 20-mer:5'-ttcaatgcaaacagtaatga-3')
            """
            ).strip(),
        )

        # More nested history
        self.assertEqual(
            ligation_product.history(),
            textwrap.dedent(
                """
            ╙── product (Dseqrecord(-18))
                └─╼ LigationSource
                    ├─╼ c (Dseqrecord(-7))
                    │   └─╼ Source
                    │       └─╼ a (Dseqrecord(-18)) ╾ Source, Source
                    ├─╼ d (Dseqrecord(-12))
                    │   └─╼ Source
                    │       └─╼  ...
                    └─╼ e (Dseqrecord(-7))
                        └─╼ Source
                            └─╼  ...
            """
            ).strip(),
        )

        self.assertEqual(
            pcr_product.history(),
            textwrap.dedent(
                f"""
            ╙── product (Dseqrecord(-22))
                └─╼ PCRSource
                    ├─╼ primer1 ({primer1.id} 8-mer:5'-ACGTACGT-3')
                    ├─╼ seq (Dseqrecord(-30))
                    └─╼ primer2 ({primer2.id} 8-mer:5'-GCGCGCGC-3')
            """
            ).strip(),
        )
        self.assertEqual(
            custom_cut_product.history(),
            textwrap.dedent(
                """
            ╙── custom_cut_product (Dseqrecord(-7))
                └─╼ Source
                    └─╼ name (Dseqrecord(-18))
            """
            ).strip(),
        )
        self.assertEqual(
            product_gateway_BP.history(),
            textwrap.dedent(
                """
            ╙── product_gateway_BP (Dseqrecord(-164))
                └─╼ GatewaySource
                    ├─╼ name (Dseqrecord(-31))
                    └─╼ name (Dseqrecord(-238))
            """
            ).strip(),
        )


class AssemblySourceTest(TestCase):

    def test_input_field_validation(self):
        source = AssemblySource(circular=True)
        self.assertEqual(source.input, [])
        self.assertRaises(
            ValidationError, AssemblySource, circular=True, input=["Hello", "World"]
        )
        self.assertRaises(ValidationError, AssemblySource, circular=True, input=None)

    def test_from_subfragment_representation(self):
        source = AssemblySource.from_subfragment_representation(
            example_assembly, example_fragments, True
        )
        self.assertEqual(
            source.input,
            [
                AssemblyFragment(
                    sequence=example_fragments[0],
                    left_location=SimpleLocation(1, 3),
                    right_location=SimpleLocation(4, 6),
                    reverse_complemented=False,
                ),
                AssemblyFragment(
                    sequence=example_fragments[1],
                    left_location=SimpleLocation(7, 9),
                    right_location=SimpleLocation(10, 12),
                    reverse_complemented=False,
                ),
            ],
        )
        self.assertEqual(source.circular, True)

    def test_to_pydantic_model(self):
        source = AssemblySource.from_subfragment_representation(
            example_assembly, example_fragments, True
        )
        self.assertEqual(
            source.to_pydantic_model(1),
            _AssemblySource(
                id=1,
                input=[
                    _AssemblyFragment(
                        sequence=id(example_fragments[0]),
                        left_location="2..3",
                        right_location="5..6",
                        reverse_complemented=False,
                    ),
                    _AssemblyFragment(
                        sequence=id(example_fragments[1]),
                        left_location="8..9",
                        right_location="11..12",
                        reverse_complemented=False,
                    ),
                ],
                circular=True,
            ),
        )


class TextFileSequenceTest(TestCase):
    def test_from_dseqrecord(self):
        dseqrecord = Dseqrecord("AAAAAAAAAAAAAAAAAAAA")
        text_file_sequence = TextFileSequence.from_dseqrecord(dseqrecord)
        self.assertEqual(
            text_file_sequence,
            TextFileSequence(
                id=id(dseqrecord),
                sequence_file_format="genbank",
                overhang_crick_3prime=0,
                overhang_watson_3prime=0,
                file_content=dseqrecord.format("genbank"),
            ),
        )

        dseqrecord = Dseqrecord(
            Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAAAAAAAA", 2, 3)
        )
        text_file_sequence = TextFileSequence.from_dseqrecord(dseqrecord)
        self.assertEqual(
            text_file_sequence,
            TextFileSequence(
                id=id(dseqrecord),
                sequence_file_format="genbank",
                overhang_crick_3prime=2,
                overhang_watson_3prime=3,
                file_content=dseqrecord.format("genbank"),
            ),
        )


class PrimerModelTest(TestCase):
    def test_from_primer(self):
        primer = Primer("AATT", name="forward_primer")
        primer_model = PrimerModel.from_primer(primer)
        self.assertEqual(
            primer_model,
            PrimerModel(id=id(primer), name="forward_primer", sequence="AATT"),
        )


class SourceInputTest(TestCase):
    def test_validation(self):
        with self.assertRaises(TypeError):
            SourceInput(sequence=1)
        with self.assertRaises(TypeError):
            SourceInput(sequence="AATT")
        with self.assertRaises(TypeError):
            SourceInput(sequence=Dseq("AA"))
        SourceInput(sequence=Dseqrecord("AA"))
        SourceInput(sequence=Primer("AATT", name="forward_primer"))


class CloningStrategyTest(TestCase):
    def test_from_dseqrecords(self):
        # If it can be serialized, it means that everything is working correctly
        # Passing the same twice does not lead to duplicates
        cs = CloningStrategy.from_dseqrecords(
            [
                golden_gate_product,
                crispr_product,
                ligation_product,
                ligation_product,
                pcr_product,
                custom_cut_product,
                product_gateway_BP,
            ]
        )
        cs.add_primer(primer1)
        cs.add_primer(guide)

        self.assertEqual(len(cs.sequences), 20)
        self.assertEqual(len(cs.primers), 3)

        # Validate that output works and that ids are reassigned
        for i in range(2):
            if i == 0:
                out_dict = cs.model_dump()
            else:
                json_str = cs.model_dump_json()
                out_dict = json.loads(json_str)
            ids_seq = {seq["id"] for seq in out_dict["sequences"]}
            ids_primers = {primer["id"] for primer in out_dict["primers"]}
            ids_sources = {source["id"] for source in out_dict["sources"]}
            all_ids = ids_seq | ids_primers | ids_sources
            self.assertEqual(max(all_ids), 23)

    def test_use_object_id_false(self):

        with id_mode(use_object_id=False):
            cs = CloningStrategy.from_dseqrecords([pcr_product])

            # Validate primer ids
            primer_ids = {primer.id for primer in cs.primers}
            self.assertEqual(primer_ids, {int(primer1.id), int(primer2.id)})

            # Validate sequence ids
            sequence_ids = {seq.id for seq in cs.sequences}
            self.assertEqual(sequence_ids, {int(pcr_product.id), int(pcr_template.id)})

            # Validate source ids
            source_ids = {source.id for source in cs.sources}
            self.assertEqual(source_ids, {int(pcr_product.id), int(pcr_template.id)})

            # Validate assembly fragment ids
            assembly_fragment_ids = {
                fragment.sequence for fragment in cs.sources[0].input
            }
            self.assertEqual(
                assembly_fragment_ids,
                {int(pcr_template.id), int(primer1.id), int(primer2.id)},
            )

            cs = CloningStrategy.from_dseqrecords([ligation_product])
            source_ids = {source.id for source in cs.sources}
            self.assertEqual(
                source_ids,
                {int(ligation_product.id), int(a.id), int(c.id), int(d.id), int(e.id)},
            )

            # Validate sequence ids
            sequence_ids = {seq.id for seq in cs.sequences}
            self.assertEqual(
                sequence_ids,
                {int(ligation_product.id), int(a.id), int(c.id), int(d.id), int(e.id)},
            )

            # Validate assembly fragment ids
            assembly_fragment_ids = {
                fragment.sequence
                for fragment in (cs.sources[0].input + cs.sources[1].input)
            }
            self.assertEqual(
                assembly_fragment_ids,
                {int(a.id), int(c.id), int(d.id), int(e.id)},
            )


class IdModeTest(TestCase):
    def test_id_mode(self):
        primer = Primer("AATT", name="forward_primer", id="123")
        primer_wrong = Primer("AATT", name="forward_primer", id="abc")
        dseqrecord = Dseqrecord("AAAAAAAAAAAAAAAAAAAA", id="456")
        dseqrecord_wrong = Dseqrecord("AAAAAAAAAAAAAAAAAAAA", id="abc")
        self.assertEqual(get_id(primer), id(primer))
        self.assertEqual(get_id(dseqrecord), id(dseqrecord))
        # The "wrong ids" are only a problem if use_object_id is False
        self.assertEqual(get_id(primer_wrong), id(primer_wrong))
        self.assertEqual(get_id(dseqrecord_wrong), id(dseqrecord_wrong))

        with id_mode(use_object_id=False):
            self.assertEqual(get_id(primer), 123)
            self.assertEqual(get_id(dseqrecord), 456)
            self.assertRaises(ValueError, get_id, primer_wrong)
            self.assertRaises(ValueError, get_id, dseqrecord_wrong)
