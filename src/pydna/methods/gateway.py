#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

from pydna.core.dseqrecord import Dseqrecord
from pydna.assembly import assembly_is_multi_site
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import GatewaySource
from Bio.SeqFeature import SimpleLocation
from pydna.methods.recombinase import Recombinase, RecombinaseCollection


def create_recombinase_dict() -> dict[str, dict[str, RecombinaseCollection]]:
    """Create a dictionary of recombinases for the Gateway reaction."""
    raw_gateway_common = {
        "attB1": "CHWVTWTgtacaaaAAANNNG",
        "attB2": "CHWVTWTgtacaagAAANNNG",
        "attB3": "CHWVTWTgtataatAAANNNG",
        "attB4": "CHWVTWTgtatagaAAANNNG",
        "attB5": "CHWVTWTgtatacaAAANNNG",
        "attL1": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtacaaaAAANNNG",
        "attL2": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtacaagAAANNNG",
        "attL3": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtataatAAANNNG",
        "attL4": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtatagaAAANNNG",
        "attL5": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtatacaAAANNNG",
        "attR1": "CHWVTWTgtacaaaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attR2": "CHWVTWTgtacaagAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attR3": "CHWVTWTgtataatAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attR4": "CHWVTWTgtatagaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attR5": "CHWVTWTgtatacaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "overlap_1": "TWTgtacaaaAAA",
        "overlap_2": "TWTgtacaagAAA",
        "overlap_3": "TWTgtataatAAA",
        "overlap_4": "TWTgtatagaAAA",
        "overlap_5": "TWTgtatacaAAA",
    }

    raw_gateway_sites_greedy = {
        **raw_gateway_common,
        "attP1": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtacaaaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attP2": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtacaagAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attP3": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtataatAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attP4": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtatagaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
        "attP5": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTgtatacaAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    }

    raw_gateway_sites_conservative = {
        **raw_gateway_common,
        "attP1": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTgtacaaaAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
        "attP2": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTgtacaagAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
        "attP3": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTgtataatAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
        "attP4": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTgtatagaAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
        "attP5": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTgtatacaAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
    }

    out_dict = {}

    for reaction in ["BP", "LR"]:
        left, right = reaction
        conservative: list[Recombinase] = []
        greedy: list[Recombinase] = []
        for i in range(1, 6):
            site1 = f"att{left}{i}"
            site2 = f"att{right}{i}"
            seq1_conservative = raw_gateway_sites_conservative[site1]
            seq2_conservative = raw_gateway_sites_conservative[site2]
            seq1_greedy = raw_gateway_sites_greedy[site1]
            seq2_greedy = raw_gateway_sites_greedy[site2]
            conservative.append(
                Recombinase(seq1_conservative, seq2_conservative, site1, site2)
            )
            greedy.append(Recombinase(seq1_greedy, seq2_greedy, site1, site2))

        out_dict[reaction] = {
            "conservative": RecombinaseCollection(conservative),
            "greedy": RecombinaseCollection(greedy),
        }
    return out_dict


recombinase_dict = create_recombinase_dict()


def gateway_overlap(
    seqx: Dseqrecord, seqy: Dseqrecord, reaction: str, greedy: bool
) -> list[tuple[int, int, int]]:
    """
    Assembly Algorithm: Find gateway overlaps. If greedy is True, it uses a more greedy consensus site to find attP sites,
    which might give false positives.

    Parameters
    ----------
    seqx : Dseqrecord
        First sequence to find overlaps.
    seqy : Dseqrecord
        Second sequence to find overlaps.
    reaction : str
        Type of Gateway reaction (BP or LR).
    greedy : bool
        If True, use greedy gateway consensus sites.

    Returns
    -------
    list[tuple[int, int, int]] A list of overlaps between the two sequences.
    """
    mode = "greedy" if greedy else "conservative"
    return recombinase_dict[reaction][mode].overlap(seqx, seqy)


def find_gateway_sites(
    seq: Dseqrecord, greedy: bool
) -> dict[str, list[SimpleLocation]]:
    """Find all gateway sites in a sequence and return a dictionary with the name and positions of the sites."""

    mode = "greedy" if greedy else "conservative"
    collection = RecombinaseCollection(
        recombinase_dict["BP"][mode].recombinases
        + recombinase_dict["LR"][mode].recombinases
    )
    return collection.find(seq)


def annotate_gateway_sites(seq: Dseqrecord, greedy: bool) -> Dseqrecord:
    """Annotate gateway sites in a sequence."""
    mode = "greedy" if greedy else "conservative"
    collection = RecombinaseCollection(
        recombinase_dict["BP"][mode].recombinases
        + recombinase_dict["LR"][mode].recombinases
    )
    return collection.annotate(seq)


def _validate(inputs, params):
    reaction_type = params.get("reaction_type")
    if reaction_type not in ["BP", "LR"]:
        raise ValueError(
            f"Invalid reaction type: {reaction_type}, can only be BP or LR"
        )


def _algorithm(params: dict):
    reaction_type = params["reaction_type"]
    greedy = params.get("greedy", False)

    def algorithm_fn(x, y, _l):
        return gateway_overlap(x, y, reaction_type, greedy)

    return algorithm_fn


def _filter_assemblies(params: dict):
    return assembly_is_multi_site if params.get("multi_site_only", False) else None


def _explain_incompatible(inputs, params):
    """Report the att sites present, so a failed reaction is diagnosable."""
    greedy = params.get("greedy", False)
    sites_in_fragments = list()
    for frag in inputs:
        sites_in_fragments.append(list(find_gateway_sites(frag, greedy).keys()))
    formatted_strings = [
        f"fragment {i + 1}: {', '.join(sites)}"
        for i, sites in enumerate(sites_in_fragments)
    ]
    raise ValueError(
        f"Inputs are not compatible for {params['reaction_type']} reaction.\n\n"
        + "\n".join(formatted_strings),
    )


gateway_assembly = register(
    Method(
        name="gateway",
        doc="""Returns the products for Gateway assembly / Gateway cloning.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    reaction_type : Literal['BP', 'LR']
        Type of Gateway reaction
    greedy : bool, optional
        If True, use greedy gateway consensus sites, by default False
    circular_only : bool, optional
        If True, only return circular assemblies, by default False
    multi_site_only : bool, optional
        If True, only return products that where 2 sites recombined. Even if input sequences
        contain multiple att sites (typically 2), a product could be generated where only one
        site recombines. That's typically not what you want, so you can set this to True to
        only return products where both att sites recombined.

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules


    Examples
    --------

    Below an example with dummy Gateway sequences, composed with minimal sequences and the consensus
    att sites.

    >>> from pydna.methods import gateway_assembly
    >>> from pydna.core.dseqrecord import Dseqrecord
    >>> attB1 = "ACAACTTTGTACAAAAAAGCAGAAG"
    >>> attP1 = "AAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    >>> attR1 = "ACAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATGCAGTCACTATG"
    >>> attL1 = "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTACAAAAAAGCAGGCT"
    >>> seq1 = Dseqrecord("aaa" + attB1 + "ccc")
    >>> seq2 = Dseqrecord("aaa" + attP1 + "ccc")
    >>> seq3 = Dseqrecord("aaa" + attR1 + "ccc")
    >>> seq4 = Dseqrecord("aaa" + attL1 + "ccc")
    >>> products_BP = gateway_assembly([seq1, seq2], "BP")
    >>> products_LR = gateway_assembly([seq3, seq4], "LR")
    >>> len(products_BP)
    2
    >>> len(products_LR)
    2

    Now let's understand the ``multi_site_only`` parameter. Let's consider a case where we are swapping fragments
    between two plasmids using an LR reaction. Experimentally, we expect to obtain two plasmids, resulting from the
    swapping between the two att sites. That's what we get if we set ``multi_site_only`` to True.

    >>> attL2 = 'aaataatgattttattttgactgatagtgacctgttcgttgcaacaaattgataagcaatgctttcttataatgccaactttgtacaagaaagctg'
    >>> attR2 = 'accactttgtacaagaaagctgaacgagaaacgtaaaatgatataaatatcaatatattaaattagattttgcataaaaaacagactacataatactgtaaaacacaacatatccagtcactatg'
    >>> insert = Dseqrecord("cccccc" + attL1 + "ccc" + attL2 + "cccccc", circular=True)
    >>> backbone = Dseqrecord("ttttt" + attR1 + "aaa" + attR2, circular=True)
    >>> products = gateway_assembly([insert, backbone], "LR", multi_site_only=True)
    >>> len(products)
    2

    However, if we set ``multi_site_only`` to False, we get 4 products, which also include the intermediate products
    where the two plasmids are combined into a single one through recombination of a single att site. This is an
    intermediate of the reaction, and typically we don't want it:

    >>> products = gateway_assembly([insert, backbone], "LR", multi_site_only=False)
    >>> print([len(p) for p in products])
    [469, 237, 232, 469]
        """,
        shape=Shape.ASSEMBLY,
        algorithm_factory=_algorithm,
        source=GatewaySource,
        limit=None,
        validate=_validate,
        filter_assemblies_factory=_filter_assemblies,
        on_empty=_explain_incompatible,
        source_fields=lambda inputs, params: {
            "reaction_type": params["reaction_type"],
            "greedy": params.get("greedy", False),
        },
        positional=(
            "frags",
            "reaction_type",
            "greedy",
            "circular_only",
            "multi_site_only",
        ),
        summary="Gateway cloning: recombine att sites in a BP or LR reaction.",
    )
)
