#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

"""Read and write SnapGene .dna files, including their cloning history.

parse_snapgene_history() reads a .dna file into a
Dseqrecord whose source holds the cloning
history. write_snapgene writes a Dseqrecord back to a .dna file.
"""

from sgffp import SgffReader, SgffObject, SgffSegment, SgffFeature, SgffWriter
from sgffp.models.history import (
    SgffHistoryNode,
    SgffHistoryNodeContent,
    SgffHistoryTreeNode,
    SgffInputSummary,
    SgffHistoryOligo,
)
from pydna.dseq import Dseq
import re
from pydna.assembly2 import (
    gibson_assembly,
    pcr_assembly,
    restriction_ligation_assembly,
    gateway_assembly,
    in_fusion_assembly,
    ligation_assembly,
    fusion_pcr_assembly,
)
from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization
from pydna.primer import Primer
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from pydna.dseqrecord import Dseqrecord
import os
from pydna.opencloning_models import (
    Source,
    AddgeneIdSource,
    UploadedFileSource,
    LigationSource,
    RestrictionAndLigationSource,
)
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Restriction import RestrictionBatch
from pydna.parsers import parse_snapgene
from pydna.utils import flatten, cutsite_to_location
import itertools
import copy
import hashlib
import json
import warnings
from sgffp.models.notes import SgffNotes

STRAND_MAP = {"+": 1, "-": -1, ".": 0, "=": 0}
UNSUPPORTED_OPERATIONS = [
    "gcCloning",
    "taCloning",
    "topoTA",
    "topoDirectional",
    "topoBlunt",
    "destroyRestrictionFragment",
]
# Steps that don't change the sequence in a way that matters for the history.
# The parser skips over them.
PASS_THROUGH_OPERATIONS = [
    "changeStrandedness",
    "editDNAEnds",
    "changeMethylation",
    "changePhosphorylation",
    "setOrigin",  # Rotation of circular sequences
    "newFileFromSelection",  # Selection of sequences from a file
]
GIBSON_LIKE_FUNCTION_DICT = {
    "gibsonAssembly": gibson_assembly,
    "inFusionCloning": in_fusion_assembly,
    "hifiAssembly": gibson_assembly,
}


class SnapgeneHistoryParserWarning(Warning):
    pass


def _segments_to_location(
    segments: list[SgffSegment], strand_int: int, seq_len: int, circular: bool
) -> SimpleLocation | CompoundLocation:
    """Convert SnapGene feature segments to a Biopython location."""
    locations = []
    for seg in segments:
        if seg.type == "gap":
            continue
        # A segment that crosses the origin is split in two
        if circular and seg.start > seg.end:
            locations.append(SimpleLocation(seg.start, seq_len, strand=strand_int))
            locations.append(SimpleLocation(0, seg.end, strand=strand_int))
        else:
            locations.append(SimpleLocation(seg.start, seg.end, strand=strand_int))

    if len(locations) == 1:
        return locations[0]
    # Biopython lists minus-strand parts in reverse order
    if strand_int == -1:
        locations = locations[::-1]
    return CompoundLocation(locations)


def _feature_to_seqfeature(
    feature: SgffFeature, seq_len: int, circular: bool
) -> SeqFeature:
    """Convert a SnapGene feature to a Biopython SeqFeature."""
    strand_int = STRAND_MAP.get(feature.strand, 0)
    location = _segments_to_location(feature.segments, strand_int, seq_len, circular)

    # Convert qualifiers: Biopython expects lists as values
    qualifiers = {
        k: [v] if not isinstance(v, list) else v for k, v in feature.qualifiers.items()
    }
    # Use the SnapGene feature name as the label. The file's own label can
    # differ (e.g. name "35S Term", label "35S ter") and can contain
    # backslashes (e.g. rep\(pMB1)), so it is kept separately to write back.
    if feature.name:
        if "label" in qualifiers:
            qualifiers[_LABEL_QUALIFIER] = qualifiers["label"]
        qualifiers["label"] = [feature.name]
    elif "label" in qualifiers:
        qualifiers["label"] = [
            v.replace("\\", " ") if isinstance(v, str) else v
            for v in qualifiers["label"]
        ]

    # Keep the feature color
    if feature.color:
        qualifiers["ApEinfo_fwdcolor"] = [feature.color]

    # Keep SnapGene settings that Biopython has no place for, such as hidden
    # features and the names and colors of each segment of a feature. They are
    # stored as JSON text so they stay with the feature when the sequence is
    # rotated or reverse complemented.
    attributes = dict(feature.extras)
    if feature.reading_frame is not None:
        attributes["readingFrame"] = str(feature.reading_frame)
    if attributes:
        qualifiers[_FEATURE_ATTRIBUTES_QUALIFIER] = [json.dumps(attributes)]

    segments = feature.segments
    segment_attributes = [_segment_attributes(s, seq_len) for s in segments]
    if any(s.extras or s.translated for s in segments) or (
        len({s.color for s in segments}) > 1
    ):
        qualifiers[_SEGMENT_ATTRIBUTES_QUALIFIER] = [json.dumps(segment_attributes)]

    return SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)


_FEATURE_ATTRIBUTES_QUALIFIER = "snapgene_feature_attributes"
_SEGMENT_ATTRIBUTES_QUALIFIER = "snapgene_segment_attributes"
_LABEL_QUALIFIER = "snapgene_label"


def _segment_attributes(segment: SgffSegment, seq_len: int) -> dict:
    """Settings of one feature segment, plus its length for matching later."""
    attributes = dict(segment.extras)
    if segment.color:
        attributes["color"] = segment.color
    if segment.translated:
        attributes["translated"] = "1"
    if segment.type == "gap":
        attributes["type"] = "gap"
    attributes["length"] = (segment.end - segment.start) % seq_len or (
        segment.end - segment.start
    )
    return attributes


def _location_to_ranges(
    location, seq_len: int, circular: bool
) -> list[tuple[int, int]]:
    """Turn a Biopython location back into SnapGene segment ranges.

    SnapGene lists segments from left to right and stores a segment that
    crosses the origin as one range with start > end. This undoes the changes
    made by _segments_to_location.
    """
    parts = [(int(p.start), int(p.end)) for p in location.parts]
    if location.strand == -1:
        parts.reverse()
    ranges: list[tuple[int, int]] = []
    for start, end in parts:
        if circular and ranges and ranges[-1][1] == seq_len and start == 0:
            ranges[-1] = (ranges[-1][0], end)  # wraps the origin
        else:
            ranges.append((start, end))
    if circular and len(ranges) > 1 and ranges[-1][1] == seq_len and ranges[0][0] == 0:
        ranges = [(ranges[-1][0], ranges[0][1])] + ranges[1:-1]

    # Make sure the segments are listed left to right. On a circular
    # sequence, segments in the right order go less than once around.
    def turn(rs):
        return sum((b[0] - a[0]) % seq_len for a, b in zip(rs, rs[1:]))

    if len(ranges) > 1 and (
        turn(ranges[::-1]) < turn(ranges) if circular else ranges[0][0] > ranges[-1][0]
    ):
        ranges.reverse()
    return ranges


def _dseq_from_seq_properties(sequence: str, circular: bool, seq_props: dict) -> Dseq:
    if circular:
        return Dseq(sequence, circular=True)
    elif (
        seq_props is not None
        and "UpstreamStickiness" in seq_props
        and "DownstreamStickiness" in seq_props
    ):
        left_ovhg = -int(seq_props.get("UpstreamStickiness"))
        right_ovhg = -int(seq_props.get("DownstreamStickiness"))
        try:
            return Dseq.from_full_sequence_and_overhangs(
                sequence, left_ovhg, right_ovhg
            )
        except ValueError as e:
            raise NotImplementedError(f"Sequence not supported: {sequence}") from e
    else:  # pragma: no cover (I don't expect this to happen)
        return Dseq(sequence)


def _history_node_to_dseqrecord(sgff_object: SgffObject, node_id: str) -> Dseqrecord:
    """Make a Dseqrecord from one step of the SnapGene history."""
    node: SgffHistoryNode = sgff_object.history.nodes[node_id]
    tree_node = sgff_object.history.get_tree_node(node_id)

    circular = tree_node.circular if tree_node else False
    seq_props = node.properties.get("AdditionalSequenceProperties")
    seq = _dseq_from_seq_properties(node.sequence, circular, seq_props)
    seq_len = node.length
    name = tree_node.name.removesuffix(".dna") if tree_node else f"node_{node_id}"
    # Spaces are kept so names match the original files

    features = []
    if seq_len != 0:
        # Steps without a sequence get one later from _source_from_tree_node
        features = [
            _feature_to_seqfeature(feat, seq_len, circular) for feat in node.features
        ]

    annotations = {}
    annotations["topology"] = "circular" if circular else "linear"
    annotations["molecule_type"] = "DNA"

    record = Dseqrecord(
        record=seq,
        id=name,
        name=name,
        description=name,
        features=features,
        annotations=annotations,
    )
    # Keep SnapGene settings for this step (such as a custom display name),
    # its notes and its primers so write_snapgene can write them back. They
    # are set after creating the record because Dseqrecord() would turn them
    # into strings.
    if tree_node and tree_node.extras:
        record.annotations["snapgene_tree_extras"] = dict(tree_node.extras)
    if node.content and node.notes.data:
        record.annotations["snapgene_node_notes"] = dict(node.notes.data)
    snapshot_primers = node.content.blocks.get(5) if node.content else None
    if snapshot_primers:
        record.annotations["snapgene_primers"] = _primers_annotation(
            snapshot_primers[0], node.sequence
        )
    if tree_node and tree_node.primers:
        record.annotations["snapgene_tree_primers"] = _primers_annotation(
            {"Primers": tree_node.primers}, node.sequence
        )
    return record


def _sequence_checksum(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def _primers_annotation(block: dict, sequence: str) -> dict:
    """Store SnapGene primers together with a checksum of their sequence.

    The primer binding sites are only valid for that exact sequence.
    """
    return {"block": copy.deepcopy(block), "checksum": _sequence_checksum(sequence)}


def _notes_to_write(record: Dseqrecord) -> dict:
    """The SnapGene notes to write for ``record``.

    These are the notes from the original file. For sequences from Addgene,
    the Addgene link the parser looks for is added if it's missing.
    """
    notes = record.annotations.get("snapgene_node_notes") or {}
    notes = dict(notes) if isinstance(notes, dict) else {}
    source = record.source
    if isinstance(source, AddgeneIdSource) and source.repository_id:
        comments = notes.get("Comments") or ""
        if not re.search(r"https://www.addgene.org/(\d+)", comments):
            link = f"https://www.addgene.org/{source.repository_id}/"
            notes["Comments"] = f"{comments} {link}".strip()
    return notes


# Primer settings that current SnapGene versions write, with their default values
_HYBRIDIZATION_PARAMS_DEFAULTS = {
    "minimumFivePrimeAnnealing": "15",
    "showAdditionalFivePrimeMatches": "1",
}


def _primers_block_for(annotation, sequence: str) -> dict | None:
    """Get the stored primers to write with sequence.

    If the sequence has changed (e.g. it was rotated), the binding sites are
    removed and SnapGene works them out again when it opens the file.
    """
    if not isinstance(annotation, dict) or not annotation.get("block"):
        return None
    block = copy.deepcopy(annotation["block"])
    # Files from old SnapGene versions lack some primer settings, and SnapGene
    # then says the file uses an old file format. Add them with SnapGene's
    # default values, keeping any values the file has.
    params = (block.get("Primers") or {}).get("HybridizationParams")
    if isinstance(params, dict):
        for key, value in _HYBRIDIZATION_PARAMS_DEFAULTS.items():
            params.setdefault(key, value)
    if annotation.get("checksum") != _sequence_checksum(sequence):
        primers = (block.get("Primers") or {}).get("Primer") or []
        for primer in primers if isinstance(primers, list) else [primers]:
            primer.pop("BindingSite", None)
    return block


def _get_restriction_batch_from_enzyme_names(
    enzyme_names: list[str],
) -> RestrictionBatch:
    if all(enz_name in rest_dict.keys() for enz_name in enzyme_names):
        return RestrictionBatch(first=[e for e in enzyme_names])
    else:  # pragma: no cover (I don't expect this to happen)
        raise ValueError(f"Unknown enzymes: {enzyme_names}")


def _get_enzyme_batch_from_input_summaries(
    input_summaries: list[SgffInputSummary],
) -> RestrictionBatch:
    enzyme_names = set(
        flatten([input_summary.enzyme_names for input_summary in input_summaries])
    )
    # Sometimes enzymes come with < > around, we remove them
    enzyme_names = set(
        enz_name.replace("<", "").replace(">", "") for enz_name in enzyme_names
    )
    # Sometimes they have Start or End tags, we remove them
    enzyme_names = enzyme_names.difference({"Start", "End"})
    return _get_restriction_batch_from_enzyme_names(enzyme_names)


def _get_sequence_inputs(source: Source) -> list[Dseqrecord]:
    """Return the starting sequences used as inputs.

    Usually these are the direct inputs. When a restriction-ligation was
    split into a digest and a ligation, the sequences before the digest are
    returned instead.
    """

    out_value = list()
    for input_value in source.input:
        if not isinstance(input_value.sequence, Dseqrecord):
            continue
        if (
            input_value.sequence.source is None
            or len(input_value.sequence.source.input) == 0
        ):
            out_value.append(input_value.sequence)
        else:
            out_value.extend(_get_sequence_inputs(input_value.sequence.source))
    return out_value


def _get_restriction_input_combinations(
    input_sequences: list[Dseqrecord], node: SgffHistoryTreeNode
) -> list[list[Dseqrecord]]:
    """All combinations of the inputs after digesting each one, to try to
    find the expected product."""

    digestion_products = list()
    for input_sequence, input_summary in zip(input_sequences, node.input_summaries):
        rb = _get_enzyme_batch_from_input_summaries([input_summary])
        digestion = input_sequence.cut(rb)
        if len(digestion) == 0:
            digestion_products.append([input_sequence])
        else:
            digestion_products.append(digestion)
    return list(itertools.product(*digestion_products))


def _source_from_tree_node(  # noqa: C901
    expected_product: Dseqrecord, node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> tuple[Source | None | int, list[SgffHistoryNode]]:
    # SnapGene doesn't always keep the sequences of a step's inputs (it marks
    # those steps as not resurrectable). Without them the step can't be
    # repeated, unless it is one that is skipped over anyway.
    if any(child.id not in sgff_object.history.nodes for child in node.children):
        if node.operation in PASS_THROUGH_OPERATIONS:
            return -1, None
        warnings.warn(
            f"Stopped at {node.operation} operation: "
            "the input sequences are not saved in the file",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []
    input_sequences = [
        _history_node_to_dseqrecord(sgff_object, child.id) for child in node.children
    ]

    expected_seguid = expected_product.seq.seguid()
    expected_dseq_and_rc = (
        expected_product.seq,
        expected_product.seq.reverse_complement(),
    )

    def parse_oligos(oligos: list[SgffHistoryOligo]) -> list[Primer]:
        return [
            Primer(oligo.sequence, name=oligo.name or f"oligo_{i + 1}")
            for i, oligo in enumerate(oligos)
        ]

    def find_expected_product(products: list[Dseqrecord]) -> Dseqrecord | None:
        if expected_product.circular:
            return next(
                (p for p in products if p.seq.seguid() == expected_seguid), None
            )
        else:
            return next((p for p in products if p.seq in expected_dseq_and_rc), None)

    if node.operation in UNSUPPORTED_OPERATIONS:
        raise NotImplementedError(f"Operation {node.operation} not supported")
    elif node.operation == "amplifyFragment":
        primers = parse_oligos(node.oligos)
        products = pcr_assembly(input_sequences[0], *primers, limit=12)
    elif node.operation == "primerDirectedMutagenesis":
        fwd_primer, *_ = parse_oligos(node.oligos)
        rvs_primer = Primer(
            fwd_primer.seq.reverse_complement(), name=f"rvs_{fwd_primer.name}"
        )
        pcr_products = pcr_assembly(
            input_sequences[0], fwd_primer, rvs_primer, limit=10
        )
        # SnapGene includes the fusion PCR in the same step
        products = list()
        for pcr_product in pcr_products:
            pcr_product.name = "mutagenesis_pcr_product"
            products.extend(fusion_pcr_assembly([pcr_product], limit=6))
    elif node.operation in PASS_THROUGH_OPERATIONS:
        return -1, None
    elif node.operation == "changeTopology":
        if expected_product.circular:
            input_sequences = [expected_product[: len(expected_product)]]
            products = ligation_assembly(input_sequences, True)
            if len(products) == 0:  # pragma: no cover (I don't expect this to happen)
                warnings.warn(
                    "Stopped at change topology operation",
                    category=SnapgeneHistoryParserWarning,
                )
                return None, []
        else:
            warnings.warn(
                "Stopped at change topology operation",
                category=SnapgeneHistoryParserWarning,
            )
            return None, []
    elif node.operation in [
        "insertFragment",
        "goldenGateAssembly",
        "insertFragments",
        "ligateFragments",
    ]:
        # Try restriction-ligation, then plain ligation, then digesting each
        # input on its own before ligating (needed when ligation would leave
        # overhangs), and finally blunt ligation.
        # When the product is circular, linear products are not made. With
        # many inputs (e.g. Golden Gate) there can be too many linear ones.
        circular_only = expected_product.circular
        rb = _get_enzyme_batch_from_input_summaries(node.input_summaries)
        if len(rb) == 0:
            products = ligation_assembly(
                input_sequences, allow_blunt=True, circular_only=circular_only
            )
        else:
            products = restriction_ligation_assembly(
                input_sequences, rb, circular_only=circular_only
            )
            if find_expected_product(products) is None:
                products = ligation_assembly(
                    input_sequences, circular_only=circular_only
                )
            if find_expected_product(products) is None:
                for combination in _get_restriction_input_combinations(
                    input_sequences, node
                ):
                    products = ligation_assembly(
                        combination, circular_only=circular_only
                    )
                    if find_expected_product(products) is not None:
                        break
                else:
                    products = ligation_assembly(
                        input_sequences, allow_blunt=True, circular_only=circular_only
                    )
    elif node.operation == "linearize":
        # This step has no input sequence, it is made from the product
        input_sequences = [expected_product.looped()]
        input_sequences[0].source = None
        rb = _get_enzyme_batch_from_input_summaries(node.input_summaries)
        if len(rb) == 0:
            warnings.warn(
                "Stopped at linearize operation without enzymes",
                category=SnapgeneHistoryParserWarning,
            )
            return None, []
        products = input_sequences[0].cut(rb)
    elif node.operation == "circularize":
        # This step has no input sequence
        assert len(input_sequences[0].seq) == 0
        history_node = sgff_object.history.nodes[node.children[0].id]
        seq_props = history_node.content.properties.get("AdditionalSequenceProperties")
        if not (
            seq_props.get("UpstreamEnzymeName")
            and seq_props.get("DownstreamEnzymeName")
        ):  # pragma: no cover (I don't expect this to happen)
            warnings.warn(
                "Stopped at circularize operation without enzymes",
                category=SnapgeneHistoryParserWarning,
            )
            return None, []
        rb = _get_restriction_batch_from_enzyme_names(
            [
                seq_props.get("UpstreamEnzymeName"),
                seq_props.get("DownstreamEnzymeName"),
            ]
        )
        original_fragments = expected_product.cut(rb)
        if (
            len(original_fragments) != 1
        ):  # pragma: no cover (I don't expect this to happen)
            warnings.warn(
                "Stopped at circularize operation not coming from a single fragment",
                category=SnapgeneHistoryParserWarning,
            )
            return None, []
        input_sequences = [original_fragments[0]]
        input_sequences[0].source = None
        products = ligation_assembly(input_sequences, allow_blunt=True)
    elif node.operation == "removeRestrictionFragment":
        rb = _get_enzyme_batch_from_input_summaries(node.input_summaries)
        products = restriction_ligation_assembly(input_sequences, rb)
        if find_expected_product(products) is None:
            raise NotImplementedError(f"Blunting not supported for {node.operation}")
    elif node.operation == "gatewayLRCloning":
        products = gateway_assembly(input_sequences, "LR")
    elif node.operation == "gatewayBPCloning":
        products = gateway_assembly(input_sequences, "BP")
    elif node.operation in ["gibsonAssembly", "inFusionCloning", "hifiAssembly"]:
        # Some inputs may be digested before the assembly
        for combination in _get_restriction_input_combinations(input_sequences, node):
            products = GIBSON_LIKE_FUNCTION_DICT[node.operation](combination, 10)
            if find_expected_product(products) is not None:
                break
    elif node.operation == "overlapFragments":
        products = fusion_pcr_assembly(input_sequences, limit=10)
    elif node.operation == "annealOligos":
        primers = parse_oligos(node.oligos)
        products = oligonucleotide_hybridization(*primers, 10)
    elif node.operation == "invalid":
        return None, []
    elif node.operation == "flip":
        # This step has no input sequence, it is made from the product
        input_sequences = [expected_product.reverse_complement()]
        input_sequences[0].source = None
        input_sequences[0].name = node.name
        products = [input_sequences[0].reverse_complement()]
    elif node.operation in [
        "remove",
        "insert",
        "replace",
        "insertReverseTranslation",
        "insertCodons",
        "insertRestrictionSite",
        "insertFeature",
    ]:
        warnings.warn(
            "Manual editing of sequences not supported",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []
    else:
        warnings.warn(
            f"Unknown operation: {node.operation}",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []

    correct_product = find_expected_product(products)

    if correct_product is None:
        raise ValueError(f"No product found for expected SEGUID {expected_seguid}")

    # Return the children nodes in the same order as the inputs
    input_index = {id(seq): i for i, seq in enumerate(input_sequences)}
    out_nodes = [node.children[input_index.get(id(seq))] for seq in input_sequences]
    return correct_product.source, out_nodes


def _is_too_many_assemblies(error: ValueError) -> bool:
    """True if the error is pydna.assembly2 refusing to list all the possible
    assemblies because there are too many."""
    return str(error).startswith("Too many assemblies")


def _replay_circular_products(
    source: RestrictionAndLigationSource | LigationSource, handle_insertion: bool
) -> list[Dseqrecord]:
    """Repeat a ligation step like source._replay_products does, but only
    making circular products."""
    input_sequences = source._get_input_sequences(handle_insertion)
    if isinstance(source, RestrictionAndLigationSource):
        return restriction_ligation_assembly(
            input_sequences, source.restriction_enzymes, circular_only=True
        )
    return ligation_assembly(
        input_sequences,
        allow_blunt=source._minimal_assembly_overlap() == 0,
        allow_partial_overlap=True,
        circular_only=True,
    )


def _can_replay_circular_only(record: Dseqrecord, error: ValueError) -> bool:
    """True if a step that failed with too many assemblies can be repeated
    making only circular products.

    Repeating a step (in normalize_history and validate_history) makes linear
    products too. For a circular product with many inputs (e.g. Golden Gate)
    there can be too many of them, although the circular ones are few.
    """
    return (
        _is_too_many_assemblies(error)
        and record.circular
        and isinstance(record.source, (RestrictionAndLigationSource, LigationSource))
    )


def _normalize_history(record: Dseqrecord) -> Dseqrecord:
    """Same as Dseqrecord.normalize_history, but ligations of a circular
    product are repeated making only circular products if there are too many
    to list."""
    if record.source is None:
        return record
    for inp in record.source.input:
        if isinstance(inp.sequence, Dseqrecord):
            inp.sequence = _normalize_history(inp.sequence)
    try:
        return record.source.normalize(record)
    except ValueError as error:
        if not _can_replay_circular_only(record, error):
            raise
    for handle_insertion in (True, False):
        products = _replay_circular_products(record.source, handle_insertion)
        try:
            product = record.source._find_product_by_seguid(record, products)
        except ValueError:
            if not handle_insertion:
                raise
            continue
        product.name = record.name
        product.id = record.id
        return product


def _validate_history(record: Dseqrecord) -> None:
    """Same as Dseqrecord.validate_history, but ligations of a circular
    product are repeated making only circular products if there are too many
    to list."""
    if record.source is None:
        return
    try:
        record.source.validate(record)
    except ValueError as error:
        if not _can_replay_circular_only(record, error):
            raise
        for handle_insertion in (True, False):
            products = _replay_circular_products(record.source, handle_insertion)
            try:
                record.source._validate_result_in_products(record, products)
                break
            except ValueError:
                if not handle_insertion:
                    raise
    for inp in record.source.input:
        if isinstance(inp.sequence, Dseqrecord):
            _validate_history(inp.sequence)


def _parse_history(
    root_record: Dseqrecord, root_node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> None:
    """Parse the history of a Dseqrecord, editing it in place."""

    repeat = True
    while repeat:
        try:
            source, out_nodes = _source_from_tree_node(
                root_record, root_node, sgff_object
            )
        except ValueError as error:
            if not _is_too_many_assemblies(error):
                raise
            warnings.warn(
                f"Stopped at {root_node.operation} operation: "
                "too many possible products to find the one in the file",
                category=SnapgeneHistoryParserWarning,
            )
            source, out_nodes = None, []
        repeat = source == -1
        if repeat:
            root_node = root_node.children[0]

    root_record.source = source
    if source is None:
        # The file's notes may say where it came from (e.g. an Addgene ID)
        if root_node.id in sgff_object.history.nodes:
            root_record.source = _source_from_metadata(
                sgff_object.history.nodes[root_node.id].content.notes
            )
        # Else, set the default source
        if root_record.source is None:
            root_record.source = _get_default_source(root_node.name)
        return
    for input_value in _get_sequence_inputs(source):
        node = out_nodes.pop(0)
        _parse_history(input_value, node, sgff_object)


def _frame_transform(
    original: Dseqrecord, target: Dseqrecord
) -> tuple[bool, int] | None:
    """How to rotate and/or reverse complement original so it reads
    exactly like target, as (reverse, shift) for _apply_frame.

    Returns None if they are not the same sequence.
    """
    if original.circular != target.circular or len(original) != len(target):
        return None
    target_str = str(target.seq).upper()
    for reverse in (False, True):
        candidate = original.reverse_complement() if reverse else original
        if not target.circular:
            if candidate.seq == target.seq:
                return reverse, 0
        else:
            shift = (str(candidate.seq).upper() * 2).find(target_str)
            if shift != -1:
                return reverse, shift
    return None


def _apply_frame(record: Dseqrecord, reverse: bool, shift: int) -> Dseqrecord:
    """Reverse complement and/or rotate record and its features."""
    if reverse:
        record = record.reverse_complement()
    if shift:
        # Dseqrecord.shifted would merge touching segments of a feature, so
        # the features are moved here instead.
        seq_len = len(record)

        def shift_location(location):
            return _segments_to_location(
                [
                    SgffSegment(
                        start=(s - shift) % seq_len,
                        end=(e - shift) % seq_len or seq_len,
                    )
                    for s, e in _location_to_ranges(location, seq_len, circular=True)
                ],
                location.strand,
                seq_len,
                circular=True,
            )

        shifted = record.shifted(shift)
        shifted.features = [
            SeqFeature(
                location=shift_location(f.location),
                type=f.type,
                id=f.id,
                qualifiers=dict(f.qualifiers),
            )
            for f in record.features
        ]
        record = shifted
    return record


def _source_from_metadata(notes: SgffNotes) -> None | Source:
    if notes.get("Comments") and (
        match := re.search(r"https://www.addgene.org/(\d+)", notes.get("Comments"))
    ):
        return AddgeneIdSource(repository_id=match.group(1))
    # This would work for sequences imported from NCBI, but GenBank files also
    # have an AccessionNumber that may be arbitrary, and there is no reliable
    # way to tell them apart without asking NCBI.
    # elif notes.get("AccessionNumber"):
    #     return NCBISequenceSource(repository_id=notes.get("AccessionNumber"))
    else:
        return None


def _get_default_source(file_name: str) -> UploadedFileSource:
    return UploadedFileSource(
        file_name=file_name,
        sequence_file_format="snapgene",
        index_in_file=0,
    )


def parse_snapgene_history(  # noqa: C901
    data: str | bytes, file_name: str = ""
) -> Dseqrecord:
    """Read a SnapGene .dna file into a Dseqrecord
    whose source holds the cloning history.

    Parameters
    ----------
    data: str | bytes
        Path to the .dna file to parse or bytes of the file content.

    file_name: str
        File name to use for the source. Taken from the path if data is a str.

    Returns
    -------
    Dseqrecord
        The sequence with its cloning history.

    Raises
    ------
    NotImplementedError
        If the file contains a cloning operation that is not yet supported.
    ValueError
        If a recorded operation cannot be reproduced (no matching product found).
    """
    root_record = parse_snapgene(data)[0]

    if isinstance(data, bytes):
        sgff_object = SgffReader.from_bytes(data)
    else:
        sgff_object = SgffReader.from_file(data)
        file_name = os.path.basename(data) if not file_name else file_name

    root_record.name = re.sub(r"\s+", "_", file_name).removesuffix(".dna")

    seq_props = sgff_object.properties.get("AdditionalSequenceProperties")
    # Biopython's reader does not handle sticky ends
    root_record.seq = _dseq_from_seq_properties(
        str(root_record.seq), root_record.circular, seq_props
    )

    # Biopython's reader changes the features (it adds primers as features,
    # uses a different label, loses colors and direction), so read them with
    # sgffp instead
    root_record.features = [
        _feature_to_seqfeature(feat, len(root_record.seq), root_record.circular)
        for feat in sgff_object.features.items
    ]

    if not sgff_object.has_history:
        root_record.source = _source_from_metadata(sgff_object.notes)
    else:
        _parse_history(root_record, sgff_object.history.tree.root, sgff_object)

    if root_record.source is None:
        root_record.source = _get_default_source(file_name)
    # _normalize_history rebuilds the record by repeating the cloning steps.
    # The result can start at a different position or be reverse complemented,
    # and its features are pieced together from the inputs. The history is only
    # valid for the rebuilt record, so it is kept, but with the file's own
    # features. write_snapgene uses "snapgene_frame" to write the sequence
    # the way the file had it.
    original_record = root_record
    root_record = _normalize_history(root_record)
    if sgff_object.has_history:
        # The records for each history step are rebuilt too, so give each one
        # the features SnapGene saved for the step with the same sequence
        snapshots: dict[str, list[Dseqrecord]] = {}
        for node_id, node in sgff_object.history.nodes.items():
            if node.length == 0:
                continue
            snapshot = _history_node_to_dseqrecord(sgff_object, node_id)
            snapshots.setdefault(snapshot.seq.seguid(), []).append(snapshot)

        def restore_features(record: Dseqrecord) -> None:
            candidates = snapshots.get(record.seq.seguid(), [])
            # Prefer a step with the same name
            for snapshot in sorted(candidates, key=lambda s: s.name != record.name):
                transform = _frame_transform(snapshot, record)
                if transform is not None:
                    record.features = _apply_frame(snapshot, *transform).features
                    for key in (
                        "snapgene_tree_extras",
                        "snapgene_node_notes",
                        "snapgene_primers",
                        "snapgene_tree_primers",
                    ):
                        if key in snapshot.annotations:
                            record.annotations.setdefault(
                                key, snapshot.annotations[key]
                            )
                    break
            if record.source is not None:
                for inp in record.source.input:
                    if isinstance(inp.sequence, Dseqrecord):
                        restore_features(inp.sequence)

        for inp in root_record.source.input:
            if isinstance(inp.sequence, Dseqrecord):
                restore_features(inp.sequence)
    transform = _frame_transform(original_record, root_record)
    if transform is not None:
        root_record.features = _apply_frame(original_record, *transform).features
        back = _frame_transform(root_record, original_record)
        if back != (False, 0):
            root_record.annotations["snapgene_frame"] = {
                "reverse": back[0],
                "shift": back[1],
                "checksum": _sequence_checksum(str(root_record.seq)),
            }
    # Keep the file's notes and primers so write_snapgene can write them back
    if sgff_object.notes.data:
        root_record.annotations["snapgene_node_notes"] = dict(sgff_object.notes.data)
    primers_block = sgff_object.blocks.get(5)
    if primers_block:
        root_record.annotations["snapgene_primers"] = _primers_annotation(
            primers_block[0], str(original_record.seq)
        )
    _validate_history(root_record)

    return root_record


_INV_STRAND_MAP = {1: "+", -1: "-", 0: ".", None: "."}
_STRAND_TO_DIRECTIONALITY = {".": "0", "+": "1", "-": "2", "=": "3"}

# SnapGene operation to write for each kind of pydna cloning step. Steps that
# SnapGene doesn't have use the closest one the parser can repeat, e.g.
# recombination is written as a Gibson assembly. Gateway cloning and a single
# sequence ligated into a circle are handled in _build_history_tree.
_SOURCE_TO_SNAPGENE_OP = {
    "PCRSource": "amplifyFragment",
    "GibsonAssemblySource": "gibsonAssembly",
    "InFusionSource": "inFusionCloning",
    "OverlapExtensionPCRLigationSource": "overlapFragments",
    "OligoHybridizationSource": "annealOligos",
    "RestrictionAndLigationSource": "insertFragments",
    "LigationSource": "ligateFragments",
    "SequenceCutSource": "linearize",
    "ReverseComplementSource": "flip",
    "AnnotationSource": "newFileFromSelection",
    "PolymeraseExtensionSource": "newFileFromSelection",
    "HomologousRecombinationSource": "gibsonAssembly",
    "CRISPRSource": "gibsonAssembly",
    "CreLoxRecombinationSource": "gibsonAssembly",
    "RecombinaseSource": "gibsonAssembly",
    "InVivoAssemblySource": "gibsonAssembly",
}

# For these cloning steps, an input that was cut with restriction enzymes is
# written as the uncut sequence plus the enzymes used, as SnapGene does.
# Otherwise the parser can't repeat ligations of sticky ends.
_COLLAPSE_CUT_PARENTS = frozenset(
    {
        "RestrictionAndLigationSource",
        "LigationSource",
        "GibsonAssemblySource",
        "InFusionSource",
        "OverlapExtensionPCRLigationSource",
        "HomologousRecombinationSource",
        "CRISPRSource",
        "InVivoAssemblySource",
        "CreLoxRecombinationSource",
        "RecombinaseSource",
    }
)


# Qualifiers that hold a feature's color. SnapGene stores the color on each
# segment instead, so these are not written as qualifiers.
_COLOR_QUALIFIER_KEYS = ("color", "ApEinfo_fwdcolor", "ApEinfo_revcolor")


def _seqfeature_to_sgff_feature(  # noqa: C901
    feature: SeqFeature, seq_len: int, circular: bool
) -> SgffFeature | None:
    """Convert a Biopython SeqFeature back to a SnapGene feature."""
    if feature.location is None:
        return None
    loc = feature.location

    segments: list[SgffSegment] = []
    for start, end in _location_to_ranges(loc, seq_len, circular):
        # A skipped region (e.g. an intron in a CDS) is a "gap" segment
        if segments and segments[-1].end < start:
            segments.append(SgffSegment(start=segments[-1].end, end=start, type="gap"))
        segments.append(SgffSegment(start=start, end=end))

    # Biopython stores qualifier values in lists, SnapGene doesn't. The label,
    # color and SnapGene settings saved by _feature_to_seqfeature are taken
    # out and used below instead of being written as qualifiers.
    qualifiers: dict = {}
    label = ""
    color = ""
    attributes: dict = {}
    segment_attributes: list[dict] = []
    for key, value in feature.qualifiers.items():
        if isinstance(value, list) and len(value) == 1:
            value = value[0]
        if key == "label":
            label = value if isinstance(value, str) else (value[0] if value else "")
        elif key in _COLOR_QUALIFIER_KEYS:
            if not color and isinstance(value, str):
                color = value
        elif key == _FEATURE_ATTRIBUTES_QUALIFIER:
            attributes = json.loads(value)
        elif key == _SEGMENT_ATTRIBUTES_QUALIFIER:
            segment_attributes = json.loads(value)
        elif key == _LABEL_QUALIFIER:
            qualifiers["label"] = value
        else:
            qualifiers[key] = value

    if color:
        for seg in segments:
            seg.color = color

    # Put back the segment settings saved by _feature_to_seqfeature. They are
    # matched to the segments by length, in reverse if the sequence was
    # reverse complemented, and skipped if the feature was changed.
    lengths = [(s.type, _segment_attributes(s, seq_len)["length"]) for s in segments]
    stored_lengths = [
        (a.get("type", "standard"), a.get("length")) for a in segment_attributes
    ]
    if stored_lengths != lengths and stored_lengths == lengths[::-1]:
        segment_attributes = segment_attributes[::-1]
        stored_lengths = lengths
    if stored_lengths == lengths:
        for seg, attrs in zip(segments, segment_attributes):
            attrs = dict(attrs)
            attrs.pop("length", None)
            attrs.pop("type", None)
            seg.color = attrs.pop("color", seg.color)
            seg.translated = attrs.pop("translated", None) == "1"
            seg.extras = attrs
    standard = [s for s in segments if s.type != "gap"]

    reading_frame = attributes.pop("readingFrame", None)

    name = label or feature.id or feature.type or ""
    if name == "<unknown id>":
        name = feature.type or ""

    return SgffFeature(
        name=name,
        type=feature.type or "misc_feature",
        strand=_INV_STRAND_MAP.get(loc.strand, "."),
        segments=segments,
        qualifiers=qualifiers,
        color=color or (standard[0].color if standard else None),
        reading_frame=int(reading_frame) if reading_frame is not None else None,
        extras=attributes,
    )


def _site_count(record: Dseqrecord, enzyme) -> int:
    """Number of cut sites enzyme has in record."""
    return len(record.seq.get_cutsites(enzyme))


def _cut_end(record: Dseqrecord, cutsite, is_left: bool) -> tuple[int, tuple[str, int]]:
    """Describe one end of a fragment cut from record.

    Returns the cut position and (enzyme name, number of sites), which is
    how SnapGene records it. An uncut end is recorded as Start/End.
    """
    if cutsite is None or cutsite[1] is None:
        if is_left:
            return 0, ("Start", 0)
        return len(record.seq), ("End", 0)
    location = cutsite_to_location(cutsite, len(record.seq))
    enzyme = cutsite[1]
    return int(location.start), (str(enzyme), _site_count(record, enzyme))


def _input_summaries_for_source(
    source,
    child_records: list[Dseqrecord],
    child_fragments: list,
    cut_edges_per_child: list[tuple | None],
) -> list[SgffInputSummary]:
    """Describe how each input was used in a cloning step.

    SnapGene records, for each input, the part of it that was used and which
    enzymes cut it. The parser needs one entry per input and the enzyme names
    to repeat the step. cut_edges_per_child holds the cut sites of inputs
    that were cut before this step (see _COLLAPSE_CUT_PARENTS).
    """
    if not child_records:
        return []
    cls = type(source).__name__

    if cls == "SequenceCutSource":
        summaries = []
        for c in child_records:
            val1, end1 = _cut_end(c, source.left_edge, True)
            val2, end2 = _cut_end(c, source.right_edge, False)
            summaries.append(
                SgffInputSummary(
                    manipulation="select", val1=val1, val2=val2, enzymes=[end1, end2]
                )
            )
        return summaries

    if cls == "PCRSource":
        return [
            SgffInputSummary(
                manipulation="amplify", val1=0, val2=len(c.seq), enzymes=[]
            )
            for c in child_records
        ]

    def fragment_end(record: Dseqrecord, location, is_left: bool):
        # Find which enzyme cut the input at this end. None if none did.
        if location is None:
            return _cut_end(record, None, is_left)
        for enzyme in source_enzymes:
            for cutsite in record.seq.get_cutsites(enzyme):
                cut_location = cutsite_to_location(cutsite, len(record.seq))
                if [(p.start, p.end) for p in cut_location.parts] == [
                    (p.start, p.end) for p in location.parts
                ]:
                    return _cut_end(record, cutsite, is_left)
        return None

    source_enzymes = []
    if cls == "RestrictionAndLigationSource":
        source_enzymes = list(getattr(source, "restriction_enzymes", []) or [])

    summaries = []
    for child, fragment, cut_edges in zip(
        child_records, child_fragments, cut_edges_per_child
    ):
        ends = None
        if cut_edges is not None:
            # The input was cut before this step
            ends = [
                _cut_end(child, cut_edges[0], True),
                _cut_end(child, cut_edges[1], False),
            ]
        elif source_enzymes:
            ends = [
                fragment_end(
                    child,
                    getattr(fragment, "left_location", None),
                    True,
                ),
                fragment_end(
                    child,
                    getattr(fragment, "right_location", None),
                    False,
                ),
            ]
            if None in ends:
                ends = None

        if ends is not None:
            (val1, end1), (val2, end2) = ends
            summaries.append(
                SgffInputSummary(
                    manipulation="insert", val1=val1, val2=val2, enzymes=[end1, end2]
                )
            )
        else:
            # The ends are unknown, so just list the enzymes
            summaries.append(
                SgffInputSummary(
                    manipulation="insert",
                    val1=0,
                    val2=len(child.seq),
                    enzymes=[(str(e), _site_count(child, e)) for e in source_enzymes],
                )
            )
    return summaries


def _make_history_node_snapshot(
    node_id: int, record: Dseqrecord, sgff_features: list[SgffFeature]
) -> SgffHistoryNode:
    """Save the sequence, features, primers and notes of one history step."""
    seq_str = str(record.seq).upper()
    seq_len = len(seq_str)
    circular = record.circular

    upstream = 0 if circular else -(record.seq.ovhg or 0)
    downstream = 0 if circular else -(record.seq.watson_ovhg or 0)
    inner_blocks: dict = {
        8: [
            {
                "AdditionalSequenceProperties": {
                    "UpstreamStickiness": str(upstream),
                    "DownstreamStickiness": str(downstream),
                    "UpstreamModification": "Unmodified",
                    "DownstreamModification": "Unmodified",
                }
            }
        ]
    }

    primers_block = _primers_block_for(
        record.annotations.get("snapgene_primers"), seq_str
    )
    if primers_block:
        inner_blocks[5] = [primers_block]

    if sgff_features:
        inner_blocks[10] = [
            {
                "features": [f.to_dict() for f in sgff_features],
                "wrapper_extras": {},
            }
        ]

    notes_dict = _notes_to_write(record)
    if notes_dict:
        inner_blocks[6] = [{"Notes": notes_dict}]

    return SgffHistoryNode(
        index=node_id,
        sequence=seq_str,
        sequence_type=1,  # compressed DNA, as SnapGene stores it
        length=seq_len,
        # The layout SnapGene uses for saved history steps
        content=SgffHistoryNodeContent({30: [inner_blocks]}),
        writer_stamp=30,
    )


def _build_history_tree(
    record: Dseqrecord,
    id_counter: list[int],
    history_nodes: dict[int, SgffHistoryNode],
    is_root: bool,
    root_name: str | None = None,
) -> SgffHistoryTreeNode:
    """Build the SnapGene history tree for record and its inputs.

    Inputs are numbered before the step that uses them, so the final sequence
    gets the highest number, as in SnapGene. The saved data for every step
    except the last is added to history_nodes.
    """
    source = record.source
    can_collapse = source is not None and type(source).__name__ in _COLLAPSE_CUT_PARENTS

    children: list[SgffHistoryTreeNode] = []
    child_records: list[Dseqrecord] = []
    child_fragments: list = []
    cut_edges_per_child: list[tuple | None] = []
    oligos: list[SgffHistoryOligo] = []

    for i, inp in enumerate(getattr(source, "input", None) or []):
        seq_obj = inp.sequence
        if isinstance(seq_obj, Primer):
            # Primers are listed on the step itself, not as inputs
            oligos.append(
                SgffHistoryOligo(
                    name=seq_obj.name or f"oligo_{i + 1}",
                    sequence=str(seq_obj.seq).upper(),
                )
            )
            continue
        if not isinstance(seq_obj, Dseqrecord):
            continue

        # Use the uncut sequence and remember the cut sites
        # (see _COLLAPSE_CUT_PARENTS)
        cut_edges = None
        cut_src = seq_obj.source
        if (
            can_collapse
            and cut_src is not None
            and type(cut_src).__name__ == "SequenceCutSource"
            and getattr(cut_src, "input", None)
        ):
            parent_seq = cut_src.input[0].sequence
            has_enzyme = any(
                edge is not None and edge[1] is not None
                for edge in (cut_src.left_edge, cut_src.right_edge)
            )
            if has_enzyme and isinstance(parent_seq, Dseqrecord):
                seq_obj = parent_seq
                cut_edges = (cut_src.left_edge, cut_src.right_edge)

        children.append(
            _build_history_tree(seq_obj, id_counter, history_nodes, is_root=False)
        )
        child_records.append(seq_obj)
        child_fragments.append(inp)
        cut_edges_per_child.append(cut_edges)

    node_id = id_counter[0]
    id_counter[0] += 1

    # Keep spaces in names so they match the original files
    name = str(root_name if is_root else (record.name or f"node_{node_id}"))

    # The features are saved twice for each step: with the step's sequence
    # and on the history tree, where SnapGene shows them.
    seq_len = len(record.seq)
    circular = bool(record.circular)
    sgff_features = [
        f
        for f in (
            _seqfeature_to_sgff_feature(feat, seq_len, circular)
            for feat in record.features
        )
        if f is not None
    ]

    def tree_feature(feat: SgffFeature) -> dict:
        # The history tree uses SnapGene's own layout for features, where
        # "Q" holds the qualifiers
        result: dict = {
            "name": feat.name,
            "type": feat.type,
            "directionality": _STRAND_TO_DIRECTIONALITY.get(feat.strand, "0"),
        }
        result.update(feat.extras)
        if feat.reading_frame is not None:
            result["readingFrame"] = str(feat.reading_frame)
        if feat.segments:
            segs = [seg.to_dict() for seg in feat.segments]
            result["Segment"] = segs[0] if len(segs) == 1 else segs
        if feat.qualifiers:

            def value(v) -> dict:
                return {"int": str(v)} if isinstance(v, int) else {"text": str(v)}

            result["Q"] = [
                {
                    "name": k,
                    "V": [value(x) for x in v] if isinstance(v, list) else value(v),
                }
                for k, v in feat.qualifiers.items()
            ]
        return result

    # The root's features are stored with the main sequence instead.
    tree_features = [] if is_root else [tree_feature(f) for f in sgff_features]

    # SnapGene settings for this step (such as a custom display name) and
    # its primers, kept by the parser
    tree_extras_val = record.annotations.get("snapgene_tree_extras") or {}
    tree_extras = dict(tree_extras_val) if isinstance(tree_extras_val, dict) else {}
    tree_primers_block = _primers_block_for(
        record.annotations.get("snapgene_tree_primers"), str(record.seq)
    )
    tree_primers = tree_primers_block["Primers"] if tree_primers_block else None

    # SnapGene operation for this step. Sequences without a cloning step get
    # "invalid". A single sequence ligated into a circle gets
    # "changeTopology", which the parser repeats with blunt ligation.
    if source is None or not getattr(source, "input", None):
        operation = "invalid"
    elif type(source).__name__ == "GatewaySource":
        reaction = str(getattr(source, "reaction_type", "LR"))
        operation = "gatewayLRCloning" if reaction == "LR" else "gatewayBPCloning"
    elif (
        type(source).__name__ == "LigationSource"
        and record.circular
        and sum(isinstance(inp.sequence, Dseqrecord) for inp in source.input) == 1
    ):
        operation = "changeTopology"
    else:
        operation = _SOURCE_TO_SNAPGENE_OP.get(type(source).__name__, "invalid")

    tree_node = SgffHistoryTreeNode(
        id=node_id,
        name=name,
        type="DNA",
        seq_len=seq_len,
        strandedness="double",
        circular=circular,
        operation=operation,
        upstream_modification="Unmodified",
        downstream_modification="Unmodified",
        resurrectable=not is_root,
        oligos=oligos,
        input_summaries=_input_summaries_for_source(
            source, child_records, child_fragments, cut_edges_per_child
        ),
        features=tree_features,
        children=children,
        primers=tree_primers,
        extras=tree_extras,
    )

    # The final sequence is saved as the main sequence of the file
    if not is_root:
        history_nodes[node_id] = _make_history_node_snapshot(
            node_id, record, sgff_features
        )

    return tree_node


def write_snapgene(record: Dseqrecord, filepath: str) -> None:
    """Write a Dseqrecord to a SnapGene .dna file.

    Writes the sequence, features and primers, and the cloning history in record.source if
    there is one.

    Each pydna cloning step is written as the closest SnapGene operation.
    Steps SnapGene doesn't have (e.g. homologous recombination) are written as
    a Gibson assembly so they can still be read back.

    The written file won't be identical to an original SnapGene file, as
    SnapGene stores some details that a Dseqrecord doesn't keep.

    Parameters
    ----------
    record:
        The Dseqrecord to write.
    filepath:
        Path of the .dna file to write.
    """
    # Write the sequence the way the original file had it, see parse_snapgene_history()
    display = record
    frame = record.annotations.get("snapgene_frame")
    # Only if the record hasn't been changed since it was read
    if isinstance(frame, dict) and frame.get("checksum") == _sequence_checksum(
        str(record.seq)
    ):
        display = _apply_frame(record, frame["reverse"], frame["shift"])

    circular = record.circular
    sequence = str(display.seq).upper()
    seq_len = len(sequence)

    sgff = SgffObject.new(
        sequence=sequence,
        topology="circular" if circular else "linear",
        strandedness="double",
        sequence_type="dna",
    )

    sgff.properties.set(
        "AdditionalSequenceProperties",
        {
            "UpstreamStickiness": str(0 if circular else -(display.seq.ovhg or 0)),
            "DownstreamStickiness": str(
                0 if circular else -(display.seq.watson_ovhg or 0)
            ),
            "UpstreamModification": "Unmodified",
            "DownstreamModification": "Unmodified",
        },
    )

    for feature in display.features:
        sgff_feature = _seqfeature_to_sgff_feature(feature, seq_len, circular)
        if sgff_feature is not None:
            sgff.features.add(sgff_feature)

    primers_block = _primers_block_for(
        record.annotations.get("snapgene_primers"), sequence
    )
    if primers_block:
        sgff.bset(5, primers_block)

    notes = _notes_to_write(record)
    if notes:
        sgff.bset(6, {"Notes": notes})

    if record.source is not None and getattr(record.source, "input", None):
        history_nodes: dict[int, SgffHistoryNode] = {}
        root_tree_node = _build_history_tree(
            record,
            id_counter=[0],
            history_nodes=history_nodes,
            is_root=True,
            root_name=os.path.basename(filepath),
        )
        sgff.bset(7, {"HistoryTree": {"Node": root_tree_node.to_dict()}})
        # Sorted so the same record always gives the same file
        sgff.bset(11, [history_nodes[k].to_dict() for k in sorted(history_nodes)])

    SgffWriter.to_file(sgff, filepath)
