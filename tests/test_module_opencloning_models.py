from pydna.opencloning_models import SequenceLocationStr
from Bio.SeqFeature import SimpleLocation
from pydna.utils import shift_location
from pydna.dseqrecord import Dseqrecord
from unittest import TestCase
from pydna.opencloning_models import (
    AssemblySource,
    AssemblyFragment,
    Source,
    SourceInput,
)
from pydna.primer import Primer
from opencloning_linkml.datamodel import (
    AssemblySource as _AssemblySource,
    AssemblyFragment as _AssemblyFragment,
    Primer as _PrimerModel,
    Source as _Source,
    SourceInput as _SourceInput,
)


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


example_assembly = [
    (1, SimpleLocation(1, 3), SimpleLocation(4, 6)),
    (2, SimpleLocation(7, 9), SimpleLocation(10, 12)),
]
example_fragments = [
    Dseqrecord("AAAAAAAAAAAAAAAAAAAA"),
    Dseqrecord("TTTTTTTTTTTTTTTTTTT"),
    Primer("AATT", name="forward_primer"),
]


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


class AssemblySourceTest(TestCase):

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
