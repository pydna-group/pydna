# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SimpleLocation
from pydna.recombinase import Recombinase


def create_recombinase_dict() -> dict[str, dict[str, list[Recombinase]]]:
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
            "conservative": conservative,
            "greedy": greedy,
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
    type = "greedy" if greedy else "conservative"
    recombinases = recombinase_dict[reaction][type]
    return sum((r.overlap(seqx, seqy) for r in recombinases), [])


def find_gateway_sites(
    seq: Dseqrecord, greedy: bool
) -> dict[str, list[SimpleLocation]]:
    """Find all gateway sites in a sequence and return a dictionary with the name and positions of the sites."""

    type = "greedy" if greedy else "conservative"
    out = dict()
    for reaction in ["BP", "LR"]:
        for rec in recombinase_dict[reaction][type]:
            out.update(rec.find(seq))
    return out


def annotate_gateway_sites(seq: Dseqrecord, greedy: bool) -> Dseqrecord:
    """Annotate gateway sites in a sequence."""
    type = "greedy" if greedy else "conservative"
    out = seq
    for reaction in ["BP", "LR"]:
        for rec in recombinase_dict[reaction][type]:
            out = rec.annotate(out)
    return seq
