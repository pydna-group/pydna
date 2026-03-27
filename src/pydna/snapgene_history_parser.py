#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Parse SnapGene .dna files and reconstruct their cloning history as
:class:`~pydna.dseqrecord.Dseqrecord` objects.

The single public entry point is :func:`parse_snapgene_history`, which reads a
``.dna`` file, resolves every step recorded in its SnapGene history, and
returns a ``Dseqrecord`` whose ``source`` attribute tree mirrors the full
cloning provenance.
"""

from sgffp import SgffReader, SgffObject, SgffSegment, SgffFeature
from sgffp.models.history import (
    SgffHistoryNode,
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
    AssemblyFragment,
    Source,
    AddgeneIdSource,
    UploadedFileSource,
)
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Restriction import RestrictionBatch
from pydna.parsers import parse_snapgene
from pydna.utils import flatten
import itertools
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
    """Convert SgffSegment list to a Biopython location."""
    locations = []
    for seg in segments:
        if seg.type == "gap":
            continue
        # Handle wrap-around for circular sequences
        if circular and seg.start > seg.end:
            locations.append(SimpleLocation(seg.start, seq_len, strand=strand_int))
            locations.append(SimpleLocation(0, seg.end, strand=strand_int))
        else:
            locations.append(SimpleLocation(seg.start, seg.end, strand=strand_int))

    if len(locations) == 1:
        return locations[0]
    # For reverse strand, reverse order
    if strand_int == -1:
        locations = locations[::-1]
    return CompoundLocation(locations)


def _feature_to_seqfeature(
    feature: SgffFeature, seq_len: int, circular: bool
) -> SeqFeature:
    """Convert an SgffFeature to a Biopython SeqFeature."""
    strand_int = STRAND_MAP.get(feature.strand, 0)
    location = _segments_to_location(feature.segments, strand_int, seq_len, circular)

    # Convert qualifiers: Biopython expects lists as values
    qualifiers = {
        k: [v] if not isinstance(v, list) else v for k, v in feature.qualifiers.items()
    }
    if feature.name and "label" not in qualifiers:
        qualifiers["label"] = [feature.name]

    return SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)


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
    """Convert a history node to a Dseqrecord.

    Sequence comes from the history node, metadata from the tree node.
    """
    node: SgffHistoryNode = sgff_object.history.nodes[node_id]
    tree_node = sgff_object.history.get_tree_node(node_id)

    circular = tree_node.circular if tree_node else False
    seq_props = node.properties.get("AdditionalSequenceProperties")
    seq = _dseq_from_seq_properties(node.sequence, circular, seq_props)
    seq_len = node.length
    name = tree_node.name if tree_node else f"node_{node_id}"
    name = re.sub(r"\s+", "_", name)  # Replace whitespace with underscores

    features = []
    if seq_len != 0:
        # Convert features from the node's content, only if sequence is present,
        # else the sequence is created by _source_from_tree_node
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
    return record


def _filter_assembly_fragments_that_are_sequences(
    input_value: list[AssemblyFragment],
) -> list[AssemblyFragment]:
    return [
        fragment
        for fragment in input_value
        if isinstance(fragment.sequence, Dseqrecord)
    ]


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
    """Return the most ancestral sequences used as inputs.

    These will typically be the immediate inputs, but when a
    restriction-ligation is split into restriction then ligation (see
    _get_restriction_input_combinations), we walk up
    the tree to find the most ancestral sequences.
    """

    out_value = list()
    for input_value in _filter_assembly_fragments_that_are_sequences(source.input):
        if (
            input_value.sequence.source is None
            or len(input_value.sequence.source.input) == 0
        ):
            out_value.append(input_value.sequence)
        else:
            out_value.extend(_get_sequence_inputs(input_value.sequence.source))
    return out_value


def _parse_oligos(oligos: list[SgffHistoryOligo]) -> list[Primer]:
    return [
        Primer(oligo.sequence, name=oligo.name or f"oligo_{i + 1}")
        for i, oligo in enumerate(oligos)
    ]


def _get_restriction_input_combinations(
    input_sequences: list[Dseqrecord], node: SgffHistoryTreeNode
) -> list[list[Dseqrecord]]:
    """Get all possible combinations of input sequences after digestion.

    This is used to try different combinations of inputs to see if we can
    find the expected product.
    """

    digestion_products = list()
    for input_sequence, input_summary in zip(input_sequences, node.input_summaries):
        rb = _get_enzyme_batch_from_input_summaries([input_summary])
        digestion = input_sequence.cut(rb)
        if len(digestion) == 0:
            digestion_products.append([input_sequence])
        else:
            digestion_products.append(digestion)
    return list(itertools.product(*digestion_products))


def _handle_circularize_operation(
    input_sequences: list[Dseqrecord],
    node: SgffHistoryTreeNode,
    sgff_object: SgffObject,
    expected_product: Dseqrecord,
) -> tuple[list[Dseqrecord], list[Dseqrecord]]:
    """Handle the circularize operation."""
    # This operation does not come with a sequence
    assert len(input_sequences[0].seq) == 0
    history_node = sgff_object.history.nodes[node.children[0].id]
    seq_props = history_node.content.properties.get("AdditionalSequenceProperties")
    if not (
        seq_props.get("UpstreamEnzymeName") and seq_props.get("DownstreamEnzymeName")
    ):  # pragma: no cover (I don't expect this to happen)
        warnings.warn(
            "Stopped at circularize operation without enzymes",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []

    rb = _get_restriction_batch_from_enzyme_names(
        [seq_props.get("UpstreamEnzymeName"), seq_props.get("DownstreamEnzymeName")]
    )
    original_fragments = expected_product.cut(rb)
    if len(original_fragments) != 1:  # pragma: no cover (I don't expect this to happen)
        warnings.warn(
            "Stopped at circularize operation not coming from a single fragment",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []

    input_sequences = [original_fragments[0]]
    input_sequences[0].source = None
    products = ligation_assembly(input_sequences, allow_blunt=True)
    return input_sequences, products


def _restriction_ligation_with_fallbacks(
    input_sequences: list[Dseqrecord],
    node: SgffHistoryTreeNode,
    find: callable,
) -> list[Dseqrecord]:
    """Try restriction-ligation, then plain ligation, then per-input digest + ligation."""
    rb = _get_enzyme_batch_from_input_summaries(node.input_summaries)
    products = restriction_ligation_assembly(input_sequences, rb)
    if find(products) is not None:
        return products
    products = ligation_assembly(input_sequences)
    if find(products) is not None:
        return products
    # Last resort: digest each input individually, then ligate
    # (happens when ligation would leave overhangs)
    for combination in _get_restriction_input_combinations(input_sequences, node):
        products = ligation_assembly(combination)
        if find(products) is not None:
            return products
    return products


def _source_from_tree_node(  # noqa: C901
    expected_product: Dseqrecord, node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> tuple[Source | None | int, list[SgffHistoryNode]]:
    input_sequences = [
        _history_node_to_dseqrecord(sgff_object, child.id) for child in node.children
    ]

    expected_seguid = expected_product.seq.seguid()
    expected_dseq_and_rc = (
        expected_product.seq,
        expected_product.seq.reverse_complement(),
    )

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
        primers = _parse_oligos(node.oligos)
        products = pcr_assembly(input_sequences[0], *primers, limit=12)
    elif node.operation == "primerDirectedMutagenesis":
        fwd_primer, *_ = _parse_oligos(node.oligos)
        rvs_primer = Primer(
            fwd_primer.seq.reverse_complement(), name=f"rvs_{fwd_primer.name}"
        )
        pcr_products = pcr_assembly(
            input_sequences[0], fwd_primer, rvs_primer, limit=10
        )
        # They bundle also the fusion pcr on the same step
        products = list()
        for pcr_product in pcr_products:
            pcr_product.name = "mutagenesis_pcr_product"
            products.extend(fusion_pcr_assembly([pcr_product], limit=6))
    elif node.operation in [
        "changeStrandedness",
        "editDNAEnds",
        "changeMethylation",
        "changePhosphorylation",
        "setOrigin",  # Rotation of circular sequences
        "newFileFromSelection",  # Selection of sequences from a file
    ]:
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
        products = _restriction_ligation_with_fallbacks(
            input_sequences, node, find_expected_product
        )
    elif node.operation == "linearize":
        # Here the input_sequence comes empty, I guess because it can be inferred
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
        input_sequences, products = _handle_circularize_operation(
            input_sequences, node, sgff_object, expected_product
        )
    elif node.operation == "removeRestrictionFragment":
        rb = _get_enzyme_batch_from_input_summaries(node.input_summaries)
        products = restriction_ligation_assembly(input_sequences, rb)
        if find_expected_product(products) is None:
            # This is using blunting as described here: https://www.neb.com/en-gb/applications/cloning-and-synthetic-biology/dna-end-modification/blunting
            raise NotImplementedError(f"Blunting not supported for {node.operation}")
    elif node.operation == "gatewayLRCloning":
        products = gateway_assembly(input_sequences, "LR")
    elif node.operation == "gatewayBPCloning":
        products = gateway_assembly(input_sequences, "BP")
    elif node.operation in ["gibsonAssembly", "inFusionCloning", "hifiAssembly"]:
        # It is possible that some inputs are digested prior to the assembly.,
        for combination in _get_restriction_input_combinations(input_sequences, node):
            products = GIBSON_LIKE_FUNCTION_DICT[node.operation](combination, 10)
            if find_expected_product(products) is not None:
                break
    elif node.operation == "overlapFragments":
        products = fusion_pcr_assembly(input_sequences, limit=10)
    elif node.operation == "annealOligos":
        primers = _parse_oligos(node.oligos)
        products = oligonucleotide_hybridization(*primers, 10)
    elif node.operation == "invalid":
        return None, []
    elif node.operation in ["remove", "insert", "replace"]:
        warnings.warn(
            "Manual editing of sequences not supported",
            category=SnapgeneHistoryParserWarning,
        )
        return None, []
    elif node.operation == "flip":
        # Here the input_sequence comes empty, I guess because it can be inferred
        input_sequences = [expected_product.reverse_complement()]
        input_sequences[0].source = None
        input_sequences[0].name = node.name
        products = [input_sequences[0].reverse_complement()]
    else:  # pragma: no cover (Tests will be updated when more are encountered)
        raise ValueError(f"Unknown operation: {node.operation}")

    correct_product = find_expected_product(products)

    if correct_product is None:
        raise ValueError(f"No product found for expected SEGUID {expected_seguid}")

    # Return the children nodes in the same order as the inputs
    input_index = {id(seq): i for i, seq in enumerate(input_sequences)}
    out_nodes = [node.children[input_index.get(id(seq))] for seq in input_sequences]
    return correct_product.source, out_nodes


def _parse_history(
    root_record: Dseqrecord, root_node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> None:
    """Parse the history of a Dseqrecord, editing it in place."""

    repeat = True
    while repeat:
        source, out_nodes = _source_from_tree_node(root_record, root_node, sgff_object)
        repeat = source == -1
        if repeat:
            root_node = root_node.children[0]

    root_record.source = source
    if source is None:
        # We may be able to get the source from the metadata (e.g. AddGene ID or NCBI accession number)
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


def _source_from_metadata(notes: SgffNotes) -> None | Source:
    if notes.get("Comments") and (
        match := re.search(r"https://www.addgene.org/(\d+)", notes.get("Comments"))
    ):
        return AddgeneIdSource(repository_id=match.group(1))
    # This works for sequences imported from NCBI in SnapGene
    # The problem is that sequences coming from genbank
    # files will also have an AccessionNumber, even if that is arbitrary (e.g. the name of
    # the sequence). I don't think there is a reliable way to distinguish between the two
    # without making a request to NCBI.
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


def parse_snapgene_history(filepath: str) -> Dseqrecord:
    """Parse a SnapGene ``.dna`` file and return a :class:`~pydna.dseqrecord.Dseqrecord`
    whose ``source`` attribute tree encodes the full cloning history.

    Parameters
    ----------
    filepath:
        Path to the ``.dna`` file to parse.

    Returns
    -------
    :class:`~pydna.dseqrecord.Dseqrecord`
        The root sequence with its provenance tree resolved.

    Raises
    ------
    NotImplementedError
        If the file contains a cloning operation that is not yet supported.
    ValueError
        If a recorded operation cannot be reproduced (no matching product found).
    """
    root_record = parse_snapgene(filepath)[0]
    sgff_object = SgffReader.from_file(filepath)

    name = os.path.basename(filepath)
    root_record.name = re.sub(r"\s+", "_", name)

    seq_props = sgff_object.properties.get("AdditionalSequenceProperties")
    # The biopython parser does not handle overhangs
    root_record.seq = _dseq_from_seq_properties(
        str(root_record.seq), root_record.circular, seq_props
    )

    if not sgff_object.has_history:
        root_record.source = _source_from_metadata(sgff_object.notes)
    else:
        _parse_history(root_record, sgff_object.history.tree.root, sgff_object)

    if root_record.source is None:
        root_record.source = _get_default_source(name)
    root_record = root_record.normalize_history()
    root_record.validate_history()
    return root_record
