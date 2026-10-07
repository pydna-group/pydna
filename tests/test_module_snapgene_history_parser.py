#!/usr/bin/env python3
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Tests for pydna.snapgene_history_parser."""

import glob
import os
import warnings

import pytest
from Bio.Restriction import BsaI
from Bio.SeqFeature import SeqFeature, SimpleLocation
from sgffp import SgffReader, SgffWriter
from pydna.dseq import Dseq
from pydna.assembly2 import restriction_ligation_assembly
from pydna.dseqrecord import Dseqrecord
from pydna.opencloning_models import (
    AddgeneIdSource,
    RestrictionAndLigationSource,
    UploadedFileSource,
)
from pydna.snapgene_history_parser import (
    parse_snapgene_history,
    write_snapgene,
    SnapgeneHistoryParserWarning,
)

TEST_FOLDER = os.path.join(os.path.dirname(__file__), "snapgene_history_files")

TEST_FILES = sorted(glob.glob(os.path.join(TEST_FOLDER, "*.dna")))

METHOD_NOT_SUPPORTED = [
    "topo_ta_cloning.dna",
    "topo_directional_cloning.dna",
    "gc_cloning_including_overhangs.dna",
    "gc_cloning_overhangs_after.dna",
    "ta_cloning_including_overhangs.dna",
    "ta_cloning_overhangs_after.dna",
    "topo_blunt_cloning.dna",
    "delete_restriction_fragment_without_compatible_overhangs.dna",
    "pcr_modifying_ends_insert_R.dna",
    "destroy_restriction_site.dna",
]

EXPECTED_VALUE_ERROR = [
    # The problem here is that it's a linear blunt ligation
    # of several sequences. pydna does not return this because
    # it would circularize.
    "blunt_linear_ligation.dna"
]


class TestSnapgeneHistoryParser:

    def test_files_exist(self):
        assert len(TEST_FILES) == 50

    @pytest.mark.parametrize("file", TEST_FILES, ids=os.path.basename)
    def test_correctly_parsed(self, file):
        filename = os.path.basename(file)
        try:
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always")
                seqr = parse_snapgene_history(file)
                file_warnings = [
                    str(w.message)
                    for w in wlist
                    if issubclass(w.category, SnapgeneHistoryParserWarning)
                ]
        except NotImplementedError as e:
            if filename not in METHOD_NOT_SUPPORTED:
                raise AssertionError(f"File {filename} not supported") from e
            return
        except ValueError as e:
            if (
                filename not in EXPECTED_VALUE_ERROR
                or "No product found for expected SEGUID" not in str(e)
            ):
                raise AssertionError(f"File {filename} not supported") from e
            return

        # Check special cases
        if filename == "import_addgene_then_clone.dna":
            assert isinstance(seqr.source.input[0].sequence.source, AddgeneIdSource)
        # Can't do this anymore, (see _source_from_metadata)
        # if filename == "import_ncbi.dna":
        #     assert isinstance(seqr.source, NCBISequenceSource)
        if filename == "import_addgene.dna":
            assert isinstance(seqr.source, AddgeneIdSource)

        # Check warnings
        if filename == "circularize_then_linearize_without_enzyme.dna":
            expected_warnings = ["Stopped at linearize operation without enzymes"]
        elif filename == "circularize.dna":
            expected_warnings = ["Stopped at change topology operation"]
        elif filename.startswith("manual_"):
            expected_warnings = ["Manual editing of sequences not supported"]
        else:
            expected_warnings = []

        assert file_warnings == expected_warnings

    def test_input_sequences_not_saved(self, tmp_path):
        # SnapGene doesn't always keep the sequences of a step's inputs. The
        # history then stops at that step with a warning instead of an error.
        sgff = SgffReader.from_file(os.path.join(TEST_FOLDER, "golden_gate.dna"))
        not_saved = sgff.history.tree.root.children[1].id
        sgff.blocks[11] = [
            node for node in sgff.blocks[11] if node["node_index"] != not_saved
        ]
        file = str(tmp_path / "golden_gate.dna")
        SgffWriter.to_file(sgff, file)

        with pytest.warns(SnapgeneHistoryParserWarning) as caught:
            seqr = parse_snapgene_history(file)
        assert [str(w.message) for w in caught] == [
            "Stopped at goldenGateAssembly operation: "
            "the input sequences are not saved in the file"
        ]
        assert isinstance(seqr.source, UploadedFileSource)

    def test_parse_snapgene_history_from_bytes(self):
        example_file = os.path.join(TEST_FOLDER, "circularize.dna")
        with open(example_file, "rb") as f:
            bytes_data = f.read()
        seqr = parse_snapgene_history(bytes_data, file_name="circularize.dna")
        assert seqr.name == "circularize"

        # Can overwrite file_name, even for files
        seqr = parse_snapgene_history(example_file, file_name="overwrite.dna")
        assert seqr.name == "overwrite"

        # Otherwise takes from file name
        seqr = parse_snapgene_history(example_file)
        assert seqr.name == "circularize"


# Files that parse_snapgene_history can read, used for the write_snapgene tests
READABLE_FILES = [
    f
    for f in TEST_FILES
    if os.path.basename(f) not in METHOD_NOT_SUPPORTED + EXPECTED_VALUE_ERROR
]

# write_snapgene can't yet write SnapGene's "circularize" step (a sticky-ended
# fragment closed into a circle), so the history of these files is not
# checked after a round trip
HISTORY_NOT_WRITTEN = ["circularize_only.dna", "rotate_restrict_rotate.dna"]


def _parse(file):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SnapgeneHistoryParserWarning)
        return parse_snapgene_history(file)


def _feature_details(features, with_positions=True):
    """Everything SnapGene stores about the features, for comparing files."""
    return sorted(
        (
            f.name,
            f.type,
            f.strand,
            f.color,
            f.reading_frame,
            sorted(f.extras.items()),
            sorted((k, str(v)) for k, v in f.qualifiers.items()),
            [
                (
                    (s.start, s.end) if with_positions else None,
                    s.type,
                    s.color,
                    s.translated,
                    sorted(s.extras.items()),
                )
                for s in f.segments
            ],
        )
        for f in features
    )


def _primer_details(primers):
    return sorted(
        (
            p.name,
            p.sequence.upper(),
            sorted((b.start, b.end, b.bound_strand) for b in p.binding_sites),
        )
        for p in primers
    )


def _history_steps(file):
    """Features of each saved history step, by step name. Positions are left
    out because steps can be written from a different start position."""
    sgff = SgffReader.from_file(file)
    steps = {}
    if sgff.has_history:
        for node_id, node in sgff.history.nodes.items():
            tree_node = sgff.history.get_tree_node(node_id)
            if tree_node and node.content and node.length:
                name = tree_node.name.removesuffix(".dna")
                steps.setdefault(name, []).append(
                    _feature_details(node.features, with_positions=False)
                )
    return steps


def _history(record):
    """The cloning steps and their sequences, as nested tuples."""
    source = record.source
    inputs = [
        inp.sequence
        for inp in (source.input if source else [])
        if isinstance(inp.sequence, Dseqrecord)
    ]
    return (
        record.seq.seguid(),
        type(source).__name__,
        tuple(_history(i) for i in inputs),
    )


class TestWriteSnapgene:

    @pytest.mark.parametrize("file", READABLE_FILES, ids=os.path.basename)
    def test_round_trip(self, file, tmp_path):
        out = str(tmp_path / "out.dna")
        record = _parse(file)
        write_snapgene(record, out)
        original = SgffReader.from_file(file)
        written = SgffReader.from_file(out)

        # Same sequence, written from the same start position
        assert written.sequence.value.upper() == original.sequence.value.upper()

        # Notes, features and primers are unchanged
        assert written.notes.data == original.notes.data
        assert _feature_details(written.features.items) == _feature_details(
            original.features.items
        )
        assert _primer_details(written.primers.items) == _primer_details(
            original.primers.items
        )

        # Each history step keeps its features
        original_steps = _history_steps(file)
        for name, features in _history_steps(out).items():
            if name in original_steps:
                assert all(f in original_steps[name] for f in features), name

    @pytest.mark.parametrize(
        "file",
        [f for f in READABLE_FILES if os.path.basename(f) not in HISTORY_NOT_WRITTEN],
        ids=os.path.basename,
    )
    def test_round_trip_keeps_history(self, file, tmp_path):
        out = str(tmp_path / "out.dna")
        record = _parse(file)
        write_snapgene(record, out)
        read_back = _parse(out)
        assert read_back.seq == record.seq
        assert _history(read_back) == _history(record)

    def test_second_round_trip_gives_same_file_contents(self, tmp_path):
        file = os.path.join(TEST_FOLDER, "gateway_single_BP_and_LR_at_once.dna")
        # Same file name, as the name is saved in the history
        (tmp_path / "first").mkdir()
        (tmp_path / "second").mkdir()
        first = str(tmp_path / "first" / "out.dna")
        second = str(tmp_path / "second" / "out.dna")
        write_snapgene(_parse(file), first)
        write_snapgene(_parse(first), second)
        assert open(first, "rb").read() == open(second, "rb").read()

    @pytest.mark.parametrize(
        "file",
        [
            "golden_gate.dna",
            "circular_ligation_3fragments.dna",
            "linear_ligation2_overhangs.dna",
        ],
    )
    def test_enzyme_cut_sites_in_history(self, file, tmp_path):
        # Both ends of each fragment are recorded with the enzyme that cut it
        # and how many sites that enzyme has, as SnapGene shows them
        file = os.path.join(TEST_FOLDER, file)
        out = str(tmp_path / "out.dna")
        write_snapgene(_parse(file), out)

        def cut_sites(path):
            step = SgffReader.from_file(path).history.tree.root
            return [(s.val1, s.val2, s.enzymes) for s in step.input_summaries]

        assert cut_sites(out) == cut_sites(file)

    @pytest.mark.parametrize("file", ["import_ncbi.dna", "import_ensembl.dna"])
    def test_hidden_features_and_segment_names(self, file, tmp_path):
        file = os.path.join(TEST_FOLDER, file)
        out = str(tmp_path / "out.dna")
        write_snapgene(_parse(file), out)

        def details(path):
            features = SgffReader.from_file(path).features.items
            hidden = sorted(f.name for f in features if f.extras.get("visible") == "0")
            segment_names = [
                [s.extras.get("name") for s in f.segments] for f in features
            ]
            return hidden, segment_names

        hidden, segment_names = details(file)
        assert hidden or any(any(names) for names in segment_names)
        assert details(out) == (hidden, segment_names)

    def test_old_primer_settings_updated(self, tmp_path):
        # Files from old SnapGene versions lack some primer settings, which
        # makes SnapGene say the file uses an old format. Missing ones are
        # added with SnapGene's defaults, the others are kept.
        record = _parse(os.path.join(TEST_FOLDER, "origin_spanning_features.dna"))
        params = record.annotations["snapgene_primers"]["block"]["Primers"][
            "HybridizationParams"
        ]
        del params["minimumFivePrimeAnnealing"]
        params["showAdditionalFivePrimeMatches"] = "0"

        out = str(tmp_path / "out.dna")
        write_snapgene(record, out)
        written = SgffReader.from_file(out).blocks[5][0]["Primers"]
        assert written["HybridizationParams"]["minimumFivePrimeAnnealing"] == "15"
        assert written["HybridizationParams"]["showAdditionalFivePrimeMatches"] == "0"

    def test_primer_binding_sites_removed_if_sequence_changes(self, tmp_path):
        file = os.path.join(TEST_FOLDER, "origin_spanning_features.dna")
        record = _parse(file)
        original = SgffReader.from_file(file).primers.items
        assert all(p.binding_sites for p in original)

        out = str(tmp_path / "out.dna")
        write_snapgene(record, out)
        assert _primer_details(SgffReader.from_file(out).primers.items) == (
            _primer_details(original)
        )

        # Rotating the record moves the binding sites, so they are left out
        # for SnapGene to find again
        rotated = record.shifted(10)
        write_snapgene(rotated, out)
        written = SgffReader.from_file(out)
        assert written.sequence.value.upper() == str(rotated.seq).upper()
        assert [(p.name, p.sequence) for p in written.primers.items] == [
            (p.name, p.sequence) for p in original
        ]
        assert not any(p.binding_sites for p in written.primers.items)

    def test_addgene_source_without_notes(self, tmp_path):
        # The Addgene link is added to the notes so the source can be read back
        record = Dseqrecord("ACGTACGTAAACCCGGGTTT", circular=True)
        record.source = AddgeneIdSource(repository_id="12345")
        out = str(tmp_path / "addgene.dna")
        write_snapgene(record, out)
        read_back = parse_snapgene_history(out)
        assert isinstance(read_back.source, AddgeneIdSource)
        assert read_back.source.repository_id == "12345"

    def test_record_without_history(self, tmp_path):
        # Linear sequence with a 4 nt 5' overhang at each end
        seq = Dseq.from_full_sequence_and_overhangs("AATTCGGATCCAAAGGTACCTTAA", 4, 4)
        record = Dseqrecord(seq, name="plain")
        record.features = [
            SeqFeature(
                SimpleLocation(4, 10, 1) + SimpleLocation(14, 20, 1),
                type="CDS",
                qualifiers={"label": ["two parts"]},
            )
        ]
        out = str(tmp_path / "plain.dna")
        write_snapgene(record, out)

        assert not SgffReader.from_file(out).has_history
        read_back = parse_snapgene_history(out)
        assert read_back.seq == record.seq
        assert isinstance(read_back.source, UploadedFileSource)
        assert read_back.source.file_name == "plain.dna"
        (feature,) = read_back.features
        assert feature.qualifiers["label"] == ["two parts"]
        assert [(p.start, p.end, p.strand) for p in feature.location.parts] == [
            (4, 10, 1),
            (14, 20, 1),
        ]

    def test_golden_gate_with_many_parts(self, tmp_path):
        # Four parts and a destination vector cut with BsaI give few circular
        # products, but more linear ones than pydna will list. Only circular
        # products are made when the product is circular.
        overhangs = ["AATG", "GCTT", "CGAA", "ACTA", "TTCG"]
        parts = [
            "acgtgactgatcgatcgtacgatgc",
            "ttgacgtagctagcatcgatgcagt",
            "gcatgcatcgtagctacgtcagtca",
            "tcagtcgatgcatgcagtcgtacgt",
        ]
        backbone = "ccatgcaagtcaggtacgtattcgagcttacgatcgatcg"
        donors = [
            Dseqrecord(
                "GGTCTCa"
                + overhangs[i]
                + part
                + overhangs[i + 1]
                + "tGAGACC"
                + backbone,
                circular=True,
                name=f"part{i + 1}",
            )
            for i, part in enumerate(parts)
        ]
        destination = Dseqrecord(
            "ttcgcatgca"
            + overhangs[0]
            + "tGAGACCaaaaaaaaaaGGTCTCa"
            + overhangs[4]
            + "ggcatcgatcgtacgtagct",
            circular=True,
            name="destination",
        )
        inputs = [destination, *donors]
        with pytest.raises(ValueError, match="Too many assemblies"):
            restriction_ligation_assembly(inputs, [BsaI])
        product = restriction_ligation_assembly(inputs, [BsaI], circular_only=True)[0]

        out = str(tmp_path / "golden_gate.dna")
        write_snapgene(product, out)
        read_back = parse_snapgene_history(out)
        assert read_back.seq.seguid() == product.seq.seguid()
        assert isinstance(read_back.source, RestrictionAndLigationSource)
        assert len(read_back.source.input) == 5

    def test_too_many_assemblies_stops_history(self, monkeypatch):
        # If a step has too many possible products to find the one in the
        # file, the history stops there instead of the file not being read
        import pydna.snapgene_history_parser as parser

        def too_many(*args, **kwargs):
            raise ValueError("Too many assemblies (99 pre-validation) to assemble")

        monkeypatch.setattr(parser, "restriction_ligation_assembly", too_many)
        file = os.path.join(TEST_FOLDER, "golden_gate.dna")
        with pytest.warns(SnapgeneHistoryParserWarning, match="too many possible"):
            record = parse_snapgene_history(file)
        assert isinstance(record.source, UploadedFileSource)
