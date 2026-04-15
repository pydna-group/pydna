# -*- coding: utf-8 -*-
import copy
from typing import Callable
from unittest import TestCase

from Bio.Restriction import BsaI, EcoRI, SalI
from Bio.Seq import reverse_complement
from Bio.SeqFeature import SimpleLocation
from opencloning_linkml.datamodel import (
    AnnotationReport,
    AssemblyFragment as _AssemblyFragment,
    AssemblySource as _AssemblySource,
    Source as _Source,
    SourceInput as _SourceInput,
    TemplateSequence,
)
from pydantic import BaseModel, ValidationError

from pydna.assembly2 import (
    cre_lox_excision_or_inversion,
    cre_lox_integration,
    crispr_integration,
    fusion_pcr_assembly,
    gateway_assembly,
    gibson_assembly,
    golden_gate_assembly,
    homologous_recombination_excision_or_inversion,
    homologous_recombination_integration,
    in_fusion_assembly,
    in_vivo_assembly,
    ligation_assembly,
    pcr_assembly,
    recombinase_excision_or_inversion,
    recombinase_integration,
)
from pydna.cre_lox import LOXP_SEQUENCE
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization
from pydna.opencloning_models import (
    AssemblyFragment,
    AssemblySource,
    CloningStrategy,
    GenomeCoordinatesSource,
    NCBISequenceSource,
    PCRSource,
    PrimerModel,
    RepositoryIdSource,
    SequenceCutSource,
    ReverseComplementSource,
    SequenceLocationStr,
    Source,
    SourceInput,
    TextFileSequence,
    UploadedFileSource,
    AnnotationSource,
    PolymeraseExtensionSource,
    _TARGET_MODEL_REGISTRY,
    _source_input_from_model,
    get_id,
    id_mode,
    read_dseqrecord_from_text_file_sequence,
)
from pydna.primer import Primer
from pydna.recombinase import Recombinase
from pydna.utils import shift_location

import glob
import json
import os
import textwrap

# Examples that will be used in several tests ==============================================
test_folder = os.path.join(os.path.dirname(__file__))


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

## Oligo hybridization

fwd_primer = Primer("ATGGC", name="fwd_primer")
rvs_primer = Primer("GCCAT", name="rvs_primer")
product_oligo_hybridization, *rest = oligonucleotide_hybridization(
    fwd_primer, rvs_primer, 3
)

product_oligo_hybridization.name = "product_oligo_hybridization"


# Recombinase example
site1 = "ATGCCCTAAaaTT"
site2 = "AAaaTTTTTTTCCCT"
recombinase = Recombinase(site1, site2)
genome = Dseqrecord(f"cccccc{site1.upper()}aaaa{site2.upper()}cccccc", name="genome")
recombinase_products = recombinase_excision_or_inversion(genome, recombinase)
recombinase_products[0].name = "excised_plasmid"
recombinase_products[1].name = "remaining_genome"
recombinase_product = recombinase_integration(
    recombinase_products[1],
    [recombinase_products[0]],
    recombinase.get_reverse_recombinase(),
)[0]
recombinase_product.name = "reconstituted_locus"


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

    def test_get_ncbi_format_coordinates(self):
        self.assertEqual(
            SequenceLocationStr("1..3").get_ncbi_format_coordinates(), (1, 3, 1)
        )
        self.assertEqual(
            SequenceLocationStr("complement(1..3)").get_ncbi_format_coordinates(),
            (1, 3, -1),
        )


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
                    │   └─╼ RestrictionEnzymeDigestionSource
                    │       └─╼ a (Dseqrecord(-18)) ╾ RestrictionEnzymeDigestionSource, RestrictionEnzymeDigestionSource
                    ├─╼ d (Dseqrecord(-12))
                    │   └─╼ RestrictionEnzymeDigestionSource
                    │       └─╼  ...
                    └─╼ e (Dseqrecord(-7))
                        └─╼ RestrictionEnzymeDigestionSource
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
                └─╼ SequenceCutSource
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
        self.assertEqual(
            product_oligo_hybridization.history(),
            textwrap.dedent(
                """
            ╙── product_oligo_hybridization (Dseqrecord(-5))
                └─╼ OligoHybridizationSource
                    ├─╼ fwd_primer (id 5-mer:5'-ATGGC-3')
                    └─╼ rvs_primer (id 5-mer:5'-GCCAT-3')
            """
            ).strip(),
        )


class AssemblySourceTest(TestCase):

    def test_input_field_validation(self):
        source = RepositoryIdSource(repository_id="1234567890")
        self.assertEqual(source.input, [])
        self.assertRaises(ValidationError, AssemblySource, circular=True, input=[])
        self.assertRaises(
            ValidationError, AssemblySource, circular=True, input=["Hello", "World"]
        )
        self.assertRaises(ValidationError, AssemblySource, circular=True, input=None)
        self.assertRaises(
            ValidationError,
            SequenceCutSource,
            left_edge=None,
            right_edge=None,
            input=[
                SourceInput(sequence=Dseqrecord("AATT")),
                SourceInput(sequence=Dseqrecord("AATT")),
            ],
        )
        self.assertRaises(
            ValidationError,
            SequenceCutSource,
            left_edge=None,
            right_edge=None,
            input=[],
        )

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

    def test_use_python_internal_id_false(self):

        with id_mode(use_python_internal_id=False):
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

            cs_dict1 = cs.model_dump()
            cs_dict2 = json.loads(cs.model_dump_json())
            for cs_dict in [cs_dict1, cs_dict2]:
                seq_ids_dict = [seq["id"] for seq in cs_dict["sequences"]]
                seq_ids_obj = [seq.id for seq in cs.sequences]
                self.assertEqual(seq_ids_dict, seq_ids_obj)


class IdModeTest(TestCase):
    def test_id_mode(self):
        primer = Primer("AATT", name="forward_primer", id="123")
        primer_wrong = Primer("AATT", name="forward_primer", id="abc")
        dseqrecord = Dseqrecord("AAAAAAAAAAAAAAAAAAAA", id="456")
        dseqrecord_wrong = Dseqrecord("AAAAAAAAAAAAAAAAAAAA", id="abc")
        self.assertEqual(get_id(primer), id(primer))
        self.assertEqual(get_id(dseqrecord), id(dseqrecord))
        # The "wrong ids" are only a problem if use_python_internal_id is False
        self.assertEqual(get_id(primer_wrong), id(primer_wrong))
        self.assertEqual(get_id(dseqrecord_wrong), id(dseqrecord_wrong))

        with id_mode(use_python_internal_id=False):
            self.assertEqual(get_id(primer), 123)
            self.assertEqual(get_id(dseqrecord), 456)
            self.assertRaises(ValueError, get_id, primer_wrong)
            self.assertRaises(ValueError, get_id, dseqrecord_wrong)


class GenomeCoordinatesSourceTest(TestCase):
    def test_coordinates_required(self):
        with self.assertRaises(ValidationError):
            GenomeCoordinatesSource(
                coordinates=None,
                repository_id="1234567890",
                assembly_accession="1234567890",
                locus_tag="1234567890",
                gene_id=1234567890,
            )
        GenomeCoordinatesSource(
            coordinates=SimpleLocation(1, 10),
            repository_id="1234567890",
            assembly_accession="1234567890",
            locus_tag="1234567890",
            gene_id=1234567890,
        )

    def test_serialize_coordinates(self):
        source = GenomeCoordinatesSource(
            coordinates=SimpleLocation(0, 10),
            repository_id="1234567890",
            assembly_accession="1234567890",
            locus_tag="1234567890",
            gene_id=1234567890,
        )
        self.assertEqual(source.model_dump()["coordinates"], "1..10")


class NCBISequenceSourceTest(TestCase):
    def test_coordinates_not_required(self):
        NCBISequenceSource(
            coordinates=None,
            repository_id="1234567890",
        )
        NCBISequenceSource(
            coordinates=SimpleLocation(1, 10),
            repository_id="1234567890",
        )


class RoundTripTest(TestCase):
    """Test that Dseqrecords -> CloningStrategy JSON -> CloningStrategy -> Dseqrecords
    produces objects with matching sequences and correct source types."""

    def _round_trip(self, products):
        cs = CloningStrategy.from_dseqrecords(products)
        json_str = cs.model_dump_json()
        cs2 = CloningStrategy.model_validate_json(json_str)
        return cs2.to_dseqrecords()

    def _assert_dseqrecord_equal(self, original: Dseqrecord, restored: Dseqrecord):
        """Assert that two Dseqrecord objects are equal."""
        restored_copy = copy.deepcopy(restored)
        restored_copy.id = original.id
        self.assertEqual(original.seq, restored_copy.seq)
        self.assertEqual(original.name, restored_copy.name)
        self.assertEqual(original.format("genbank"), restored_copy.format("genbank"))

    def _assert_source_type_matches(self, original, restored):
        """Assert the source type is the same (or both None)."""
        if original.source is None:
            self.assertIsNone(restored.source)
        else:
            self.assertIsNotNone(restored.source)
            self.assertEqual(type(original.source), type(restored.source))

    def test_golden_gate_round_trip(self):
        results = self._round_trip([golden_gate_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(golden_gate_product, results[0])
        self._assert_source_type_matches(golden_gate_product, results[0])

    def test_crispr_round_trip(self):
        results = self._round_trip([crispr_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(crispr_product, results[0])
        self._assert_source_type_matches(crispr_product, results[0])

    def test_ligation_round_trip(self):
        results = self._round_trip([ligation_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(ligation_product, results[0])
        self._assert_source_type_matches(ligation_product, results[0])

    def test_pcr_round_trip(self):
        results = self._round_trip([pcr_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(pcr_product, results[0])
        self._assert_source_type_matches(pcr_product, results[0])

    def test_custom_cut_round_trip(self):
        results = self._round_trip([custom_cut_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(custom_cut_product, results[0])
        self._assert_source_type_matches(custom_cut_product, results[0])

    def test_gateway_round_trip(self):
        results = self._round_trip([product_gateway_BP])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(product_gateway_BP, results[0])
        self._assert_source_type_matches(product_gateway_BP, results[0])

    def test_oligo_hybridization_round_trip(self):
        results = self._round_trip([product_oligo_hybridization])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(product_oligo_hybridization, results[0])
        self._assert_source_type_matches(product_oligo_hybridization, results[0])

    def test_recombinase_round_trip(self):
        results = self._round_trip([recombinase_product])
        self.assertEqual(len(results), 1)
        self._assert_dseqrecord_equal(recombinase_product, results[0])
        self._assert_source_type_matches(recombinase_product, results[0])

    def test_multiple_products_round_trip(self):
        products = [
            golden_gate_product,
            crispr_product,
            ligation_product,
            pcr_product,
            custom_cut_product,
            product_gateway_BP,
        ]
        results = self._round_trip(products)
        self.assertEqual(len(results), len(products))
        for orig, restored in zip(products, results):
            self._assert_dseqrecord_equal(orig, restored)
            self._assert_source_type_matches(orig, restored)

    def test_registry_covers_all_source_subclasses(self):
        """Verify the registry has entries for all Source subclasses with ClassVar TARGET_MODEL."""

        def _collect_subclasses(cls):
            result = set()
            for sub in cls.__subclasses__():
                result.add(sub)
                result.update(_collect_subclasses(sub))
            return result

        all_subclasses = _collect_subclasses(Source)
        for cls in all_subclasses:
            target = cls.__dict__.get("TARGET_MODEL")
            if target is not None and not isinstance(target, property):
                self.assertIn(
                    target,
                    _TARGET_MODEL_REGISTRY,
                    f"{cls.__name__}'s TARGET_MODEL ({target}) not in registry",
                )


class ValidateTest(TestCase):
    """Test Source.validate() and Dseqrecord.validate_history()."""

    def test_validate_golden_gate(self):
        golden_gate_product.validate_history()

    def test_validate_crispr(self):
        crispr_product.validate_history()

    def test_validate_ligation(self):
        ligation_product.validate_history()

    def test_validate_pcr(self):
        pcr_product.validate_history()

    def test_validate_cut(self):
        custom_cut_product.validate_history()

    def test_validate_gateway(self):
        product_gateway_BP.validate_history()

    def test_validate_oligo_hybridization(self):
        product_oligo_hybridization.validate_history()

    def test_validate_recombinase(self):
        recombinase_product.validate_history()

    def test_validate_custom_cut(self):
        custom_cut_product.validate_history()

    def test_validate_non_recursive(self):
        copy_ligation_product = copy.deepcopy(ligation_product)
        first_fragment = copy_ligation_product.source.input[0].sequence
        first_fragment.source.input[0].sequence = Dseqrecord("ATGC")
        with self.assertRaises(ValueError):
            copy_ligation_product.validate_history()
        copy_ligation_product.validate_history(recursive=False)

    def test_validate_reverse_complement(self):
        seq = Dseqrecord("ATGCATGC")
        rc = seq.reverse_complement()
        self.assertIsInstance(rc.source, ReverseComplementSource)
        rc.source.validate(rc)

    def test_validate_restriction_enzyme_cut(self):
        """Validate a restriction enzyme cut (SequenceCutSource with enzyme)."""
        template = Dseqrecord("aaGAATTCcc")
        products = template.cut(EcoRI)
        for p in products:
            p.validate_history()

    def test_validate_no_source(self):
        """validate_history on a Dseqrecord with no source should return silently."""
        seq = Dseqrecord("ATGC")
        seq.validate_history()

    def test_validate_no_inputs(self):
        """Source with no inputs (external) should return silently."""
        source = UploadedFileSource(
            file_name="test.gb", index_in_file=0, sequence_file_format="genbank"
        )
        source.validate(Dseqrecord("ATGC"))

    def test_validate_wrong_sequence_raises(self):
        """Mutating the sequence should cause validate to raise ValueError."""

        template = Dseqrecord("aaGAATTCcc")
        products = template.cut(EcoRI)
        product = products[0]
        # Create a wrong result with a different sequence
        wrong = Dseqrecord("aaTTTAA")
        wrong.source = product.source
        with self.assertRaises(ValueError):
            wrong.source.validate(wrong)

        # Also works "high up" in the history
        copy_golden_gate = copy.deepcopy(golden_gate_product)
        copy_golden_gate.source.input[0].sequence = Dseqrecord("aaTTTAA")
        with self.assertRaises(ValueError):
            copy_golden_gate.validate_history()

        # Works with primers as well
        copy_pcr = copy.deepcopy(pcr_product)
        copy_pcr.source.input[0].sequence = Primer("AATT")
        with self.assertRaises(ValueError):
            copy_pcr.validate_history()

        # Same for CRISPR
        copy_crispr = copy.deepcopy(crispr_product)
        for i, inp in enumerate(copy_crispr.source.input):
            if isinstance(inp.sequence, Primer):
                index_of_primer = i
                break

        copy_crispr.source.input[index_of_primer].sequence = Primer("AATT")
        with self.assertRaises(ValueError):
            copy_crispr.validate_history()

        # Same for hybridization
        copy_hybridization = copy.deepcopy(product_oligo_hybridization)
        copy_hybridization.source.input[0].sequence = Primer("AATT")
        with self.assertRaises(ValueError):
            copy_hybridization.validate_history()

    def test_validate_annotation(self):
        source = AnnotationSource(
            annotation_tool="plannotate",
            input=[SourceInput(sequence=Dseqrecord("ATGC"))],
        )
        annotated_seq = Dseqrecord("ATGC")
        annotated_seq.add_feature(0, 3, type_="gene", label="gene1")
        source.validate(annotated_seq)

        # error for different sequence
        source2 = AnnotationSource(
            annotation_tool="plannotate",
            input=[SourceInput(sequence=Dseqrecord("ATG"))],
        )
        with self.assertRaises(ValueError) as e:
            source2.validate(annotated_seq)

        self.assertIn(
            "AnnotationSource input sequence does not match result", e.exception.args[0]
        )

        # error for primer
        source3 = AnnotationSource(
            annotation_tool="plannotate",
            input=[SourceInput(sequence=Primer("AATT"))],
        )
        with self.assertRaises(ValueError) as e:
            source3.validate(annotated_seq)
        self.assertIn(
            "AnnotationSource input must be a Dseqrecord", e.exception.args[0]
        )

    def test_validate_polymerase_extension(self):
        prev_seq = Dseqrecord(Dseq.from_full_sequence_and_overhangs("ATGC", -1, -1))
        new_seq = Dseqrecord("ATGC")
        new_seq.source = PolymeraseExtensionSource(
            input=[SourceInput(sequence=prev_seq)]
        )
        new_seq.validate_history()

        # error for primer
        new_seq.source.input[0].sequence = Primer("AATT")
        with self.assertRaises(ValueError) as e:
            new_seq.validate_history()
        self.assertIn(
            "PolymeraseExtensionSource input must be a Dseqrecord", e.exception.args[0]
        )

    def test_validate_homologous_recombination(self):
        homology = "AAGTCCGTTCGTTTTACCTG"
        genome = Dseqrecord(f"aaaaaa{homology}cccc", name="genome")
        insert_seq = Dseqrecord(f"{homology}tttt{homology}", name="insert")
        products = homologous_recombination_integration(genome, [insert_seq], 20)
        products = homologous_recombination_excision_or_inversion(products[0], 20)
        for p in products:
            p.validate_history()

    def test_validate_cre_lox_excision(self):
        genome = Dseqrecord(
            f"cccccc{LOXP_SEQUENCE}aaaa{LOXP_SEQUENCE}cccccc", name="genome"
        )
        products = cre_lox_excision_or_inversion(genome)
        for p in products:
            p.validate_history()

    def test_validate_cre_lox_integration(self):
        linear = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
        circular = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)
        products = cre_lox_integration(linear, [circular])
        for p in products:
            p.validate_history()

    def test_validate_gibson_like(self):
        homology1 = "GAGTCTCC"
        homology2 = "TCAGAAGT"
        homology3 = "TTCTTCAG"
        fragments = [
            Dseqrecord(f"{homology1}acgatAAtgctcc{homology2}", name="f1"),
            Dseqrecord(f"{homology2}tcatGGGG{homology3}", name="f2"),
            Dseqrecord(f"{homology3}atataTTTT{homology1}", name="f3"),
        ]
        for func in [
            gibson_assembly,
            in_fusion_assembly,
            fusion_pcr_assembly,
            in_vivo_assembly,
        ]:
            products = func(fragments, limit=8)
            for p in products:
                p.validate_history()

    def test_validate_examples_opencloning(self):

        for file in glob.glob(f"{test_folder}/examples_opencloning/*.json"):
            with open(file, "r") as f:
                data = json.load(f)
            cloning_strategy = CloningStrategy.model_validate(data)
            for product in cloning_strategy.to_dseqrecords():
                product.validate_history()

    def test_validate_cloning_strategy(self):
        cs = CloningStrategy.from_dseqrecords([ligation_product])
        cs.validate()
        copy_ligation_product = copy.deepcopy(ligation_product)
        for inp in copy_ligation_product.source.input:
            inp.sequence = inp.sequence.reverse_complement()
        cs2 = CloningStrategy.from_dseqrecords([copy_ligation_product])
        with self.assertRaises(ValueError):
            cs2.validate()


def _replace_sequence(
    cloning_strategy: CloningStrategy,
    sequence_id: int,
    modifier: Callable[[Dseqrecord], Dseqrecord],
) -> CloningStrategy:
    """Utility function for tests"""
    model = next((s for s in cloning_strategy.sequences if s.id == sequence_id))
    seqr = read_dseqrecord_from_text_file_sequence(model)
    new_seqr = modifier(seqr)
    new_model = TextFileSequence.from_dseqrecord(new_seqr)
    new_model.id = model.id
    cs_out = copy.deepcopy(cloning_strategy)
    cs_out.sequences.remove(model)
    cs_out.sequences.append(new_model)
    return cs_out


class NormalizeTest(TestCase):

    def _common_testing_function(self, product: Dseqrecord):
        product_seguid = product.seq.seguid()
        cs = CloningStrategy.from_dseqrecords([product])

        def modify_seq(seq: Dseqrecord) -> Dseqrecord:
            if seq.circular:
                return seq.shifted(15).reverse_complement()
            else:
                return seq.reverse_complement()

        ids2replace = cs.get_ids_of_sequences_that_are_inputs()
        for seq_id in ids2replace:
            cs = _replace_sequence(cs, seq_id, modify_seq)
        seqr_with_wrong_history = cs.to_dseqrecords()
        self.assertEqual(len(seqr_with_wrong_history), 1)
        with self.assertRaises(ValueError):
            seqr_with_wrong_history[0].validate_history()

        normalized_seqr = seqr_with_wrong_history[0].normalize_history()
        self.assertEqual(normalized_seqr.seq.seguid(), product_seguid)
        normalized_seqr.validate_history()

    def test_golden_gate(self):
        self._common_testing_function(golden_gate_product)

    def test_crispr(self):
        self._common_testing_function(crispr_product)

    def test_ligation(self):
        self._common_testing_function(ligation_product)

    def test_pcr(self):
        self._common_testing_function(pcr_product)

    def test_custom_cut(self):
        with self.assertRaises(ValueError):
            self._common_testing_function(custom_cut_product)

    def test_gateway(self):
        self._common_testing_function(product_gateway_BP)

    def test_polymerase_extension(self):
        prev_seq = Dseqrecord(Dseq.from_full_sequence_and_overhangs("ATTT", -1, -1))
        new_seq = Dseqrecord("ATTT")
        new_seq.source = PolymeraseExtensionSource(
            input=[SourceInput(sequence=prev_seq)]
        )
        self._common_testing_function(new_seq)

    def test_normalize_examples_opencloning(self):
        for file in glob.glob(f"{test_folder}/examples_opencloning/*.json"):
            with open(file, "r") as f:
                data = json.load(f)
            cloning_strategy = CloningStrategy.model_validate(data)
            for product in cloning_strategy.to_dseqrecords():
                self._common_testing_function(product)

    def test_normalize_error(self):
        dummy_sequence = Dseqrecord("ATGC")
        copy_golden_gate = copy.deepcopy(golden_gate_product)
        copy_golden_gate.source.input[0].sequence = dummy_sequence
        with self.assertRaises(ValueError):
            copy_golden_gate.normalize_history()

    def test_normalize_annotation_source(self):
        for circular in [True, False]:
            seq = Dseqrecord("ATTTT", circular=circular)
            annotation_source = AnnotationSource(
                annotation_tool="plannotate",
                input=[SourceInput(sequence=seq)],
                annotation_report=[AnnotationReport()],
            )
            source_no_report = copy.deepcopy(annotation_source)
            source_no_report.annotation_report = []
            mods = [
                copy.deepcopy(seq),
                seq.reverse_complement(),
            ]
            if circular:
                mods.append(seq.shifted(15).reverse_complement())
            for i, mod in enumerate(mods):
                mod.source = annotation_source
                newseq = mod.normalize_history()
                self.assertEqual(newseq.seq, seq.seq)
                if i == 0:
                    self.assertEqual(newseq.source, annotation_source)
                else:
                    self.assertEqual(newseq.source, source_no_report)
                self.assertEqual(newseq.source.input[0].sequence.seq, seq.seq)

    def test_normalize_annotation_source_errors(self):
        seq = Dseqrecord("ATTTT")
        annotation_source = AnnotationSource(
            annotation_tool="plannotate",
            input=[SourceInput(sequence=Dseqrecord("AA"))],
            annotation_report=[AnnotationReport()],
        )
        with self.assertRaises(ValueError) as e:
            annotation_source.normalize(seq)
        self.assertEqual(
            str(e.exception), "AnnotationSource input sequence does not match result"
        )

        annotation_source.input[0].sequence = Primer("AATT")
        with self.assertRaises(ValueError) as e:
            annotation_source.normalize(seq)
        self.assertIn("AnnotationSource input must be a Dseqrecord", str(e.exception))

    def test_normalize_cloning_strategy(self):
        # Just tests that this function just calls normalize_history on the sequences
        # It uses ligation_product because for all file_contents to match, the ids must
        # be assigned from the beginning.
        with id_mode(use_python_internal_id=False):
            cs = CloningStrategy.from_dseqrecords([ligation_product])
            copy_ligation_product = copy.deepcopy(ligation_product)
            for inp in copy_ligation_product.source.input:
                this_id = inp.sequence.id
                inp.sequence = inp.sequence.reverse_complement()
                inp.sequence.source = None
                inp.sequence.id = this_id
            cs_wrong = CloningStrategy.from_dseqrecords([copy_ligation_product])
            cs_norm = cs_wrong.normalize()
            cs_norm2 = CloningStrategy.from_dseqrecords(
                [s.normalize_history() for s in cs_wrong.to_dseqrecords()]
            )
            self.assertEqual(cs_norm, cs_norm2)
            self.assertEqual(cs, cs.normalize())
            self.assertNotEqual(cs_norm, cs)

            # Same with PCR products, to ensure that it works with primers as well
            cs = CloningStrategy.from_dseqrecords([pcr_product])
            copy_pcr_product = copy.deepcopy(pcr_product)
            inp = copy_pcr_product.source.input[1]
            this_id = inp.sequence.id
            inp.sequence = inp.sequence.reverse_complement()
            inp.sequence.source = None
            inp.sequence.id = this_id
            cs_wrong = CloningStrategy.from_dseqrecords([copy_pcr_product])
            cs_norm = cs_wrong.normalize()
            cs_norm2 = CloningStrategy.from_dseqrecords(
                [s.normalize_history() for s in cs_wrong.to_dseqrecords()]
            )
            self.assertEqual(cs_norm, cs_norm2)
            self.assertEqual(cs, cs.normalize())
            self.assertNotEqual(cs_norm, cs)


class MiscTests(TestCase):
    """Miscellaneous tests to complete coverage."""

    def test_source_input_from_model_error(self):
        with self.assertRaises(ValueError):
            _source_input_from_model(_SourceInput(sequence=1), {2: "blah"}, {3: "bluh"})

    def test_source_raises_replay_products_error(self):
        with self.assertRaises(NotImplementedError):
            Source()._replay_products()

    def test_minimal_assembly_overlap_error(self):
        # If an incomplete assembly is passed for validation.
        source = PCRSource(
            input=[SourceInput(sequence=Dseqrecord("ATGC"))], circular=False
        )
        with self.assertRaises(ValueError) as e:
            source._minimal_assembly_overlap()
        self.assertEqual(str(e.exception), "Assembly is not complete")

    def test_to_dseqrecords_fails_with_template_sequence(self):
        with id_mode(use_python_internal_id=False):
            cs = CloningStrategy.from_dseqrecords([Dseqrecord("ATGC", id="1")])
            cs.sequences.append(TemplateSequence(id=23))
            with self.assertRaises(NotImplementedError):
                cs.to_dseqrecords()

    def test_to_dseqrecords_fails_with_missing_source(self):
        with id_mode(use_python_internal_id=False):
            cs = CloningStrategy.from_dseqrecords([Dseqrecord("ATGC", id="1")])
            cs.sources[0].id = 2
            with self.assertRaises(ValueError) as e:
                cs.to_dseqrecords()
            self.assertEqual(str(e.exception), "Missing source for sequence 1")

    def test_replay_products_error(self):
        source = AssemblySource(
            input=[SourceInput(sequence=Dseqrecord("ATGC"))], circular=False
        )
        with self.assertRaises(NotImplementedError):
            source._replay_products()
