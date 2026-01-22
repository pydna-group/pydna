# -*- coding: utf-8 -*-
"""
Improved implementation of the assembly module. To see a list of issues with the previous implementation,
see [issues tagged with fixed-with-new-assembly-model](https://github.com/pydna-group/pydna/issues?q=is%3Aissue%20state%3Aopen%20label%3Afixed-with-new-assembly-model)
"""

import networkx as nx
import itertools
from Bio.SeqFeature import SimpleLocation, Location

from Bio.Restriction.Restriction import RestrictionBatch
import regex
import copy

from pydna.utils import (
    shift_location,
    flatten,
    location_boundaries,
    locations_overlap,
    sum_is_sticky,
    limit_iterator,
    create_location,
)
from pydna._pretty import pretty_str as ps
from pydna.common_sub_strings import common_sub_strings as common_sub_strings_str
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.primer import Primer
from pydna.seqrecord import SeqRecord
from pydna.types import (
    CutSiteType,
    # TODO: allow user to enforce multi-site
    EdgeRepresentationAssembly,
    SubFragmentRepresentationAssembly,
    AssemblyAlgorithmType,
    SequenceOverlap,
    AssemblyEdgeType,
)
from pydna.gateway import gateway_overlap, find_gateway_sites
from pydna.cre_lox import cre_loxP_overlap
from pydna.alphabet import anneal_strands

from typing import TYPE_CHECKING, Callable, Literal
from pydna.opencloning_models import (
    AssemblySource,
    RestrictionAndLigationSource,
    GibsonAssemblySource,
    InFusionSource,
    OverlapExtensionPCRLigationSource,
    InVivoAssemblySource,
    LigationSource,
    GatewaySource,
    HomologousRecombinationSource,
    CreLoxRecombinationSource,
    PCRSource,
    SourceInput,
    CRISPRSource,
)
from pydna.crispr import cas9
import warnings

if TYPE_CHECKING:  # pragma: no cover
    from Bio.Restriction import AbstractCut


def gather_overlapping_locations(
    locs: list[Location], fragment_length: int
) -> list[tuple[Location, ...]]:
    """
    Turn a list of locations into a list of tuples of those locations, where each tuple contains
    locations that overlap. For example, if locs = [loc1, loc2, loc3], and loc1 and loc2 overlap,
    the output will be [(loc1, loc2), (loc3,)].
    """
    # Make a graph with all the locations as nodes
    G = nx.Graph()
    for i, loc in enumerate(locs):
        G.add_node(i, location=loc)

    # Add edges between nodes that overlap
    for i in range(len(locs)):
        for j in range(i + 1, len(locs)):
            if locations_overlap(locs[i], locs[j], fragment_length):
                G.add_edge(i, j)

    # Get groups of overlapping locations
    groups = list()
    for loc_set in nx.connected_components(G):
        groups.append(tuple(locs[i] for i in loc_set))

    # Sort by location of the first element in each group (does not matter which since they are overlapping)
    groups.sort(key=lambda x: location_boundaries(x[0])[0])

    return groups


def ends_from_cutsite(
    cutsite: CutSiteType, seq: Dseq
) -> tuple[tuple[str, str], tuple[str, str]]:
    """Get the sticky or blunt ends created by a restriction enzyme cut.

    Parameters
    ----------
    cutsite : CutSiteType
        A tuple ((cut_watson, ovhg), enzyme) describing where the cut occurs
    seq : _Dseq
        The DNA sequence being cut

    Raises
    ------
    ValueError
        If cutsite is None

    Returns
    -------
    tuple[tuple[str, str], tuple[str, str]]
        A tuple of two tuples, each containing the type of end ('5\'', '3\'', or 'blunt')
        and the sequence of the overhang. The first tuple is for the left end, second for the right end.

    >>> from Bio.Restriction import NotI
    >>> x = Dseq("ctcgGCGGCCGCcagcggccg")
    >>> x.get_cutsites(NotI)
    [((6, -4), NotI)]
    >>> ends_from_cutsite(x.get_cutsites(NotI)[0], x)
    (("5'", 'ggcc'), ("5'", 'ggcc'))
    """

    if cutsite is None:
        raise ValueError("None is not supported")

    cut_watson, cut_crick, ovhg = seq.get_cut_parameters(cutsite, is_left=None)
    if ovhg < 0:
        # TODO check the edge in circular
        return (
            ("5'", str(seq[cut_watson:cut_crick].reverse_complement()).lower()),
            ("5'", str(seq[cut_watson:cut_crick]).lower()),
        )
    elif ovhg > 0:
        return (
            ("3'", str(seq[cut_crick:cut_watson]).lower()),
            ("3'", str(seq[cut_crick:cut_watson].reverse_complement()).lower()),
        )

    return ("blunt", ""), ("blunt", "")


def restriction_ligation_overlap(
    seqx: Dseqrecord,
    seqy: Dseqrecord,
    enzymes=RestrictionBatch,
    partial=False,
    allow_blunt=False,
) -> list[SequenceOverlap]:
    """Assembly algorithm to find overlaps that would result from restriction and ligation.

    Like in sticky and gibson, the order matters (see example below of partial overlap)

    Parameters
    ----------
    seqx : Dseqrecord
        The first sequence
    seqy : Dseqrecord
        The second sequence
    enzymes : RestrictionBatch
        The enzymes to use
    partial : bool
        Whether to allow partial overlaps
    allow_blunt : bool
        Whether to allow blunt ends

    Returns
    -------
    list[SequenceOverlap]
        A list of overlaps between the two sequences

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import restriction_ligation_overlap
    >>> from Bio.Restriction import EcoRI, RgaI, DrdI, EcoRV
    >>> x = Dseqrecord("ccGAATTCaa")
    >>> y = Dseqrecord("aaaaGAATTCgg")
    >>> restriction_ligation_overlap(x, y, [EcoRI])
    [(3, 5, 4)]
    >>> restriction_ligation_overlap(y, x, [EcoRI])
    [(5, 3, 4)]

    Partial overlap, note how it is not symmetric

    >>> x = Dseqrecord("GACTAAAGGGTC")
    >>> y = Dseqrecord("AAGCGATCGCAAGCGATCGCAA")
    >>> restriction_ligation_overlap(x, y, [RgaI, DrdI], partial=True)
    [(6, 5, 1), (6, 15, 1)]
    >>> restriction_ligation_overlap(y, x, [RgaI, DrdI], partial=True)
    []

    Blunt overlap, returns length of the overlap 0

    >>> x = Dseqrecord("aaGATATCcc")
    >>> y = Dseqrecord("ttttGATATCaa")
    >>> restriction_ligation_overlap(x, y, [EcoRV], allow_blunt=True)
    [(5, 7, 0)]
    >>> restriction_ligation_overlap(y, x, [EcoRV], allow_blunt=True)
    [(7, 5, 0)]

    """
    cuts_x = seqx.seq.get_cutsites(*enzymes)
    cuts_y = seqy.seq.get_cutsites(*enzymes)
    # If blunt ends are allowed, something similar to this could be done to allow
    # joining with linear sequence ends, but for now it messes up with the only_adjacent_edges
    # case
    # if allow_blunt:
    #     if not seqx.circular:
    #         cuts_x.append(((len(seqx), 0), None))
    #     if not seqy.circular:
    #         cuts_y.append(((0, 0), None))
    matches = list()
    for cut_x, cut_y in itertools.product(cuts_x, cuts_y):
        # A blunt end
        if allow_blunt and cut_x[0][1] == cut_y[0][1] == 0:
            matches.append((cut_x[0][0], cut_y[0][0], 0))
            continue

        # Otherwise, test overhangs
        overlap = sum_is_sticky(
            ends_from_cutsite(cut_x, seqx.seq)[0],
            ends_from_cutsite(cut_y, seqy.seq)[1],
            partial,
        )
        if not overlap:
            continue
        x_watson, x_crick, x_ovhg = seqx.seq.get_cut_parameters(cut_x, is_left=False)
        y_watson, y_crick, y_ovhg = seqy.seq.get_cut_parameters(cut_y, is_left=True)
        # Positions where the overlap would start for full overlap
        left_x = x_watson if x_ovhg < 0 else x_crick
        left_y = y_watson if y_ovhg < 0 else y_crick

        # Correct por partial overlaps
        left_x += abs(x_ovhg) - overlap

        matches.append((left_x, left_y, overlap))
    return matches


def combine_algorithms(*algorithms: AssemblyAlgorithmType) -> AssemblyAlgorithmType:
    """
    Combine assembly algorithms, if any of them returns a match, the match is returned.

    This can be used for example in a ligation where you want to allow both sticky and blunt end ligation.
    """

    def combined(seqx, seqy, limit):
        matches = list()
        for algorithm in algorithms:
            matches += algorithm(seqx, seqy, limit)
        return matches

    return combined


def blunt_overlap(
    seqx: Dseqrecord, seqy: Dseqrecord, limit=None
) -> list[SequenceOverlap]:
    """
    Assembly algorithm to find blunt overlaps. Used for blunt ligation.

    It basically returns [(len(seqx), 0, 0)] if the right end of seqx is blunt and the
    left end of seqy is blunt (compatible with blunt ligation). Otherwise, it returns an empty list.

    Parameters
    ----------
    seqx : Dseqrecord
        The first sequence
    seqy : Dseqrecord
        The second sequence
    limit : int
        There for compatibility, but it is ignored

    Returns
    -------
    list[SequenceOverlap]
        A list of overlaps between the two sequences

    >>> from pydna.assembly2 import blunt_overlap
    >>> from pydna.dseqrecord import Dseqrecord
    >>> x = Dseqrecord("AAAAAA")
    >>> y = Dseqrecord("TTTTTT")
    >>> blunt_overlap(x, y)
    [(6, 0, 0)]
    """
    if (
        seqx.seq.three_prime_end()[0] == "blunt"
        and seqy.seq.five_prime_end()[0] == "blunt"
    ):
        return [(len(seqx), 0, 0)]
    return []


def common_sub_strings(
    seqx: Dseqrecord, seqy: Dseqrecord, limit=25
) -> list[SequenceOverlap]:
    """
    Assembly algorithm to find common substrings of length == limit. see the docs of
    the function common_sub_strings_str for more details. It is case insensitive.

    >>> from pydna.dseqrecord import Dseqrecord
    >>> x = Dseqrecord("TAAAAAAT")
    >>> y = Dseqrecord("CCaAaAaACC")
    >>> common_sub_strings(x, y, limit=5)
    [(1, 2, 6), (1, 3, 5), (2, 2, 5)]
    """
    query_seqx = str(seqx.seq).upper()
    query_seqy = str(seqy.seq).upper()
    if seqx.circular:
        query_seqx = query_seqx * 2
    if seqy.circular:
        query_seqy = query_seqy * 2
    results = common_sub_strings_str(query_seqx, query_seqy, limit)

    if not seqx.circular and not seqy.circular:
        return results

    # Remove matches that start on the second copy of the sequence
    if seqx.circular:
        results = [r for r in results if r[0] < len(seqx)]
    if seqy.circular:
        results = [r for r in results if r[1] < len(seqy)]

    # Trim lengths that span more than the sequence
    if seqx.circular or seqy.circular:
        max_match_length = min(len(seqx), len(seqy))
        results = [(r[0], r[1], min(r[2], max_match_length)) for r in results]

    # Edge case where the sequences are identical
    if len(seqx.seq) == len(seqy.seq):
        full_match = next((r for r in results if r[2] == len(seqx.seq)), None)
        if full_match is not None:
            return [full_match]

    # Remove duplicate matches, see example below
    # Let's imagine the following two sequences, where either seqy or both are circular
    # seqx: 01234
    # seqy: 123450, circular
    #
    # common_sub_strings would return [(0, 5, 5), (1, 0, 4)]
    # Actually, (1, 0, 4) is a subset of (0, 5, 5), the part
    # that does not span the origin. To remove matches like this,
    # We find matches where the origin is spanned in one of the sequences
    # only, and then remove the subset of that match that does not span the origin.
    shifted_matches = set()
    for x, y, length in results:
        x_span_origin = seqx.circular and x + length > len(seqx)
        y_span_origin = seqy.circular and y + length > len(seqy)
        if x_span_origin and not y_span_origin:
            shift = len(seqx) - x
            shifted_matches.add((0, y + shift, length - shift))
        elif not x_span_origin and y_span_origin:
            shift = len(seqy) - y
            shifted_matches.add((x + shift, 0, length - shift))
    return [r for r in results if r not in shifted_matches]


def _get_trim_end_info(
    end_info: tuple[str, str], trim_ends: str, is_five_prime: bool
) -> int | None:
    """Utility function to get the trim information for terminal_overlap."""
    if end_info[0] == trim_ends:
        return len(end_info[1]) if is_five_prime else len(end_info[1]) * -1
    return 0 if is_five_prime else None


def terminal_overlap(
    seqx: Dseqrecord, seqy: Dseqrecord, limit=25, trim_ends: None | str = None
):
    """
    Assembly algorithm to find terminal overlaps (e.g. for Gibson assembly).
    The order matters, we want alignments like:

    ::

        seqx:    oooo------xxxx
        seqy:              xxxx------oooo
        Product: oooo------xxxx------oooo

        Not like:

        seqx:               oooo------xxxx
        seqy:     xxxx------oooo
        Product (unwanted): oooo

    Parameters
    ----------
    seqx : Dseqrecord
        The first sequence
    seqy : Dseqrecord
        The second sequence
    limit : int
        Minimum length of the overlap
    trim_ends : str
        The ends to trim, either '5' or '3'
        If None, no trimming is done

    Returns
    -------
    list[SequenceOverlap]
        A list of overlaps between the two sequences

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import terminal_overlap
    >>> x = Dseqrecord("ttactaAAAAAA")
    >>> y = Dseqrecord("AAAAAAcgcacg")
    >>> terminal_overlap(x, y, limit=5)
    [(6, 0, 6), (7, 0, 5)]
    >>> terminal_overlap(y, x, limit=5)
    []

    Trimming the ends:
    >>> from pydna.dseq import Dseq
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import terminal_overlap
    >>> x = Dseqrecord(Dseq.from_full_sequence_and_overhangs("aaaACGT", 0, 3))
    >>> y = Dseqrecord(Dseq.from_full_sequence_and_overhangs("ACGTccc", 3, 0))
    >>> terminal_overlap(x, y, limit=4)
    [(3, 0, 4)]
    >>> terminal_overlap(x, y, limit=4, trim_ends="5'")
    [(3, 0, 4)]
    >>> terminal_overlap(x, y, limit=4, trim_ends="3'")
    []
    """

    if trim_ends is not None and trim_ends not in ["5'", "3'"]:
        raise ValueError("trim_ends must be '5' or '3'")

    if trim_ends is None:
        trim_x_left, trim_x_right, trim_y_left, trim_y_right = (0, None, 0, None)
        stringx = str(seqx.seq).upper()
        stringy = str(seqy.seq).upper()
    else:
        trim_x_right = _get_trim_end_info(
            seqx.seq.three_prime_end(), trim_ends, is_five_prime=False
        )
        trim_y_left = _get_trim_end_info(
            seqy.seq.five_prime_end(), trim_ends, is_five_prime=True
        )

        # I actually don't think these two are needed, since only the terminal
        # join between x_right and y_left is tested, but maybe there is some edge-case
        # that I am missing, so keeping them just in case.
        trim_x_left = _get_trim_end_info(
            seqx.seq.five_prime_end(), trim_ends, is_five_prime=True
        )
        trim_y_right = _get_trim_end_info(
            seqy.seq.three_prime_end(), trim_ends, is_five_prime=False
        )

        stringx = str(seqx.seq[trim_x_left:trim_x_right]).upper()
        stringy = str(seqy.seq[trim_y_left:trim_y_right]).upper()

    # We have to convert to list because we need to modify the matches
    matches = [
        list(m)
        for m in common_sub_strings_str(stringx, stringy, limit)
        if (m[1] == 0 and m[0] + m[2] == len(stringx))
    ]

    # Shift the matches if the left end has been trimmed
    for match in matches:
        match[0] += trim_x_left
        match[1] += trim_y_left

    # convert to tuples again
    return [tuple(m) for m in matches]


def gibson_overlap(seqx: Dseqrecord, seqy: Dseqrecord, limit=25):
    """
    Assembly algorithm to find terminal overlaps for Gibson assembly.
    It is a wrapper around terminal_overlap with trim_ends="5'".
    """

    return terminal_overlap(seqx, seqy, limit, trim_ends="5'")


def in_fusion_overlap(seqx: Dseqrecord, seqy: Dseqrecord, limit=25):
    """
    Assembly algorithm to find terminal overlaps for in-fusion assembly.
    It is a wrapper around terminal_overlap with trim_ends="3'".
    """
    return terminal_overlap(seqx, seqy, limit, trim_ends="3'")


def pcr_fusion_overlap(seqx: Dseqrecord, seqy: Dseqrecord, limit=25):
    """
    Assembly algorithm to find terminal overlaps for PCR fusion assembly.
    It is a wrapper around terminal_overlap with trim_ends=None.
    """
    return terminal_overlap(seqx, seqy, limit, trim_ends=None)


def sticky_end_sub_strings(seqx: Dseqrecord, seqy: Dseqrecord, limit: bool = False):
    """
    Assembly algorithm for ligation of sticky ends.

    For now, if limit 0 / False (default) only full overlaps are considered.
    Otherwise, partial overlaps are also returned.

    Parameters
    ----------
    seqx : Dseqrecord
        The first sequence
    seqy : Dseqrecord
        The second sequence
    limit : bool
        Whether to allow partial overlaps

    Returns
    -------
    list[SequenceOverlap]
        A list of overlaps between the two sequences


    Ligation of fully overlapping sticky ends, note how the order matters

    >>> from pydna.dseq import Dseq
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import sticky_end_sub_strings
    >>> x = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAAAA", 0, 3))
    >>> y = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAAAA", 3, 0))
    >>> sticky_end_sub_strings(x, y, limit=False)
    [(3, 0, 3)]
    >>> sticky_end_sub_strings(y, x, limit=False)
    []

    Ligation of partially overlapping sticky ends, specified with limit=True

    >>> x = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAAAA", 0, 2))
    >>> y = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAAAA", 3, 0))
    >>> sticky_end_sub_strings(x, y, limit=False)
    []
    >>> sticky_end_sub_strings(x, y, limit=True)
    [(4, 0, 2)]

    """

    overlap = sum_is_sticky(
        seqx.seq.three_prime_end(), seqy.seq.five_prime_end(), limit
    )
    if overlap:
        return [(len(seqx) - overlap, 0, overlap)]
    return []


def zip_match_leftwards(
    seqx: SeqRecord, seqy: SeqRecord, match: SequenceOverlap
) -> SequenceOverlap:
    """
    Starting from the rightmost edge of the match, return a new match encompassing the max
    number of bases. This can be used to return a longer match if a primer aligns for longer
    than the limit or a shorter match if there are mismatches. This is convenient to maintain
    as many features as possible. It is used in PCR assembly.

    >>> seq = Dseqrecord('AAAAACGTCCCGT')
    >>> primer = Dseqrecord('ACGTCCCGT')
    >>> match = (13, 9, 0) # an empty match at the end of each
    >>> zip_match_leftwards(seq, primer, match)
    (4, 0, 9)

    Works in circular molecules if the match spans the origin:
    >>> seq = Dseqrecord('TCCCGTAAAAACG', circular=True)
    >>> primer = Dseqrecord('ACGTCCCGT')
    >>> match = (6, 9, 0)
    >>> zip_match_leftwards(seq, primer, match)
    (10, 0, 9)

    """

    query_x = seqrecord2_uppercase_DNA_string(seqx)
    query_y = seqrecord2_uppercase_DNA_string(seqy)

    # In circular sequences, the match may go beyond the left-most edge of the sequence if it spans
    # the origin:
    # Primer:          ACGTCCCGT
    #                  |||||||||
    # Circular seq:    ACGTCCCGT -> Equivalent to Dseqrecord('CCCGTACGT', circular=True)
    #                      ^
    #                      Origin
    # We would start from the last T and move leftwards, but we would stop at the origin
    # For those cases we shift by length, then go back

    end_on_x = match[0] + match[2]
    if isinstance(seqx, Dseqrecord) and seqx.circular and end_on_x <= len(seqx):
        end_on_x += len(seqx)

    end_on_y = match[1] + match[2]
    if isinstance(seqy, Dseqrecord) and seqy.circular and end_on_y <= len(seqy):
        end_on_y += len(seqy)

    count = 0
    for x, y in zip(reversed(query_x[:end_on_x]), reversed(query_y[:end_on_y])):
        if x != y:
            break
        count += 1

    # Shift back by length if needed
    start_on_x = (end_on_x - count) % len(seqx)
    start_on_y = (end_on_y - count) % len(seqy)

    return (start_on_x, start_on_y, count)


def zip_match_rightwards(
    seqx: Dseqrecord, seqy: Dseqrecord, match: SequenceOverlap
) -> SequenceOverlap:
    """Same as zip_match_leftwards, but towards the right."""

    query_x = seqrecord2_uppercase_DNA_string(seqx)
    query_y = seqrecord2_uppercase_DNA_string(seqy)

    start_on_x, start_on_y, _ = match
    count = 0
    for x, y in zip(query_x[start_on_x:], query_y[start_on_y:]):
        if x != y:
            break
        count += 1
    return (start_on_x, start_on_y, count)


def seqrecord2_uppercase_DNA_string(seqr: SeqRecord) -> str:
    """
    Transform a Dseqrecord to a sequence string where U is replaced by T, everything is upper case and
    circular sequences are repeated twice. This is used for PCR, to support primers with U's (e.g. for USER cloning).
    """
    out = str(seqr.seq).upper().replace("U", "T")
    if isinstance(seqr, Dseqrecord) and seqr.circular:
        return out * 2
    return out


def primer_template_overlap(
    seqx: Dseqrecord | Primer, seqy: Dseqrecord | Primer, limit=25, mismatches=0
) -> list[SequenceOverlap]:
    """
    Assembly algorithm to find overlaps between a primer and a template. It accepts mismatches.
    When there are mismatches, it only returns the common part between the primer and the template.

    If seqx is a primer and seqy is a template, it represents the binding of a forward primer.
    If seqx is a template and seqy is a primer, it represents the binding of a reverse primer,
    where the primer has been passed as its reverse complement (see examples).

    Parameters
    ----------
    seqx : Dseqrecord | Primer
        The primer
    seqy : Dseqrecord | Primer
        The template
    limit : int
        Minimum length of the overlap
    mismatches : int
        Maximum number of mismatches (only substitutions, no deletion or insertion)

    Returns
    -------
    list[SequenceOverlap]
        A list of overlaps between the primer and the template

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.primer import Primer
    >>> from pydna.assembly2 import primer_template_overlap
    >>> template = Dseqrecord("AATTAGCAGCGATCGAGT", circular=True)
    >>> primer = Primer("TTAGCAGC")
    >>> primer_template_overlap(primer, template, limit=8, mismatches=0)
    [(0, 2, 8)]

    This actually represents the binding of the primer ``GCTGCTAA`` (reverse complement)
    >>> primer_template_overlap(template, primer, limit=8, mismatches=0)
    [(2, 0, 8)]
    >>> primer_template_overlap(primer, template.reverse_complement(), limit=8, mismatches=0)
    []
    >>> primer_template_overlap(primer.reverse_complement(), template, limit=8, mismatches=0)
    []
    """

    if isinstance(seqx, Primer) and isinstance(seqy, Dseqrecord):
        primer = seqx
        template = seqy
        reverse_primer = False
    elif isinstance(seqx, Dseqrecord) and isinstance(seqy, Primer):
        primer = seqy
        template = seqx
        reverse_primer = True
    else:
        raise ValueError(
            "One of the sequences must be a primer and the other a Dseqrecord"
        )

    if len(primer) < limit:
        return []

    subject = seqrecord2_uppercase_DNA_string(template)
    query = (
        seqrecord2_uppercase_DNA_string(primer[:limit])
        if reverse_primer
        else seqrecord2_uppercase_DNA_string(primer[-limit:])
    )

    re_matches = list(
        regex.finditer(
            "(" + query + "){s<=" + str(mismatches) + "}", subject, overlapped=True
        )
    )
    re_matches += list(
        regex.finditer(
            "(?r)(" + query + "){s<=" + str(mismatches) + "}", subject, overlapped=True
        )
    )

    out = set()
    for re_match in re_matches:

        start, end = re_match.span()

        # For circular sequences the same match is returned twice unless it falls
        # on the origin, we eliminate duplicates here
        if start >= len(template):
            continue

        # This extends match beyond the limit if the primer aligns more than that
        # and reduces the match if the primer has mismatches
        if reverse_primer:
            # Match in the same format as other assembly algorithms
            starting_match = (start, 0, end - start)
            out.add(zip_match_rightwards(template, primer, starting_match))
        else:
            # Match in the same format as other assembly algorithms
            starting_match = (len(primer) - limit, start, end - start)
            out.add(zip_match_leftwards(primer, template, starting_match))

    return list(sorted(out))


def reverse_complement_assembly(
    assembly: EdgeRepresentationAssembly, fragments: list[Dseqrecord]
) -> EdgeRepresentationAssembly:
    """Complement an assembly, i.e. reverse the order of the fragments and the orientation of the overlaps."""
    new_assembly = list()
    for u, v, locu, locv in assembly:
        f_u = fragments[abs(u) - 1]
        f_v = fragments[abs(v) - 1]
        new_assembly.append((-v, -u, locv._flip(len(f_v)), locu._flip(len(f_u))))
    return new_assembly[::-1]


def filter_linear_subassemblies(
    linear_assemblies: list[EdgeRepresentationAssembly],
    circular_assemblies: list[EdgeRepresentationAssembly],
    fragments: list[Dseqrecord],
) -> list[EdgeRepresentationAssembly]:
    """Remove linear assemblies which are sub-assemblies of circular assemblies"""
    all_circular_assemblies = circular_assemblies + [
        reverse_complement_assembly(c, fragments) for c in circular_assemblies
    ]
    filtered_assemblies = [
        assem
        for assem in linear_assemblies
        if not any(is_sublist(assem, c, True) for c in all_circular_assemblies)
    ]
    # I don't think the line below is necessary, but just in case
    # filtered_assemblies = [l for l in filtered_assemblies if not any(is_sublist(reverse_complement_assembly(l, fragments), c, True) for c in all_circular_assemblies)]
    return filtered_assemblies


def remove_subassemblies(
    assemblies: list[EdgeRepresentationAssembly],
) -> list[EdgeRepresentationAssembly]:
    """Filter out subassemblies, i.e. assemblies that are contained within another assembly.

    For example:
        [(1, 2, '1[8:14]:2[1:7]'), (2, 3, '2[10:17]:3[1:8]')]
        [(1, 2, '1[8:14]:2[1:7]')]
    The second one is a subassembly of the first one.
    """

    # Sort by length, longest first
    assemblies = sorted(assemblies, key=len, reverse=True)

    filtered_assemblies = list()
    for assembly in assemblies:
        # Check if this assembly is a subassembly of any of the assemblies we have already found
        if not any(is_sublist(assembly, a) for a in filtered_assemblies):
            filtered_assemblies.append(assembly)

    return filtered_assemblies


def assembly2str(assembly: EdgeRepresentationAssembly) -> str:
    """Convert an assembly to a string representation, for example:
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    becomes:
    ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]')

    The reason for this is that by default, a feature '[8:14]' when present in a tuple
    is printed to the console as ``SimpleLocation(ExactPosition(8), ExactPosition(14), strand=1)`` (very long).
    """
    return str(tuple(f"{u}{lu}:{v}{lv}" for u, v, lu, lv in assembly))


def assembly2str_tuple(assembly: EdgeRepresentationAssembly) -> str:
    """Convert an assembly to a string representation, like
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    """
    return str(tuple((u, v, str(lu), str(lv)) for u, v, lu, lv in assembly))


def assembly_has_mismatches(
    fragments: list[Dseqrecord], assembly: EdgeRepresentationAssembly
) -> bool:
    """Check if an assembly has mismatches. This should never happen and if so it returns an error."""
    for u, v, loc_u, loc_v in assembly:
        seq_u = fragments[u - 1] if u > 0 else fragments[-u - 1].reverse_complement()
        seq_v = fragments[v - 1] if v > 0 else fragments[-v - 1].reverse_complement()
        # TODO: Check issue where extraction failed, and whether it would give problems here
        if (
            str(loc_u.extract(seq_u).seq).upper()
            != str(loc_v.extract(seq_v).seq).upper()
        ):
            return True
    return False


def assembly_is_circular(
    assembly: EdgeRepresentationAssembly, fragments: list[Dseqrecord]
) -> bool:
    """
    Based on the topology of the locations of an assembly, determine if it is circular.
    This does not work for insertion assemblies, that's why assemble takes the optional argument is_insertion.
    """
    if assembly[0][0] != assembly[-1][1]:
        return False
    elif (
        isinstance(fragments[abs(assembly[0][0]) - 1], Dseqrecord)
        and fragments[abs(assembly[0][0]) - 1].circular
    ):
        return True
    else:
        return (
            location_boundaries(assembly[0][2])[0]
            > location_boundaries(assembly[-1][3])[0]
        )


def assemble(
    fragments: list[Dseqrecord],
    assembly: EdgeRepresentationAssembly,
    is_insertion: bool = False,
) -> Dseqrecord:
    """Generate a Dseqrecord from an assembly and a list of fragments."""

    if is_insertion:
        is_circular = False
    else:
        is_circular = assembly_is_circular(assembly, fragments)

    subfragment_representation = edge_representation2subfragment_representation(
        assembly, is_circular
    )

    # Sanity check
    for asm_edge in assembly:
        u, v, loc_u, loc_v = asm_edge
        f_u = fragments[u - 1] if u > 0 else fragments[-u - 1].reverse_complement()
        f_v = fragments[v - 1] if v > 0 else fragments[-v - 1].reverse_complement()
        seq_u = str(loc_u.extract(f_u).seq)
        seq_v = str(loc_v.extract(f_v).seq.rc())
        # Test if seq_u and seq_v anneal
        if not anneal_strands(seq_u, seq_v):
            raise ValueError("Mismatch in assembly")

    # We transform into Dseqrecords (for primers)
    dseqr_fragments = [
        f if isinstance(f, Dseqrecord) else Dseqrecord(f) for f in fragments
    ]
    subfragments = get_assembly_subfragments(
        dseqr_fragments, subfragment_representation
    )

    # Length of the overlaps between consecutive assembly fragments
    fragment_overlaps = [len(e[-1]) for e in assembly]
    out_dseqrecord = subfragments.pop(0)

    for fragment, overlap in zip(subfragments, fragment_overlaps):
        out_dseqrecord.seq = out_dseqrecord.seq.cast_to_ds_right()
        out_dseqrecord.seq = out_dseqrecord.seq.exo1_end(overlap)
        fragment.seq = fragment.seq.cast_to_ds_left()
        fragment.seq = fragment.seq.exo1_front(overlap)
        out_dseqrecord += fragment

    # For circular assemblies, process the fragment and loop
    if is_circular:
        out_dseqrecord.seq = out_dseqrecord.seq.cast_to_ds_left()
        out_dseqrecord.seq = out_dseqrecord.seq.cast_to_ds_right()
        overlap = fragment_overlaps[-1]
        out_dseqrecord.seq = out_dseqrecord.seq.exo1_front(overlap)
        out_dseqrecord.seq = out_dseqrecord.seq.exo1_end(overlap)
        out_dseqrecord = out_dseqrecord.looped()

    out_dseqrecord.source = AssemblySource.from_subfragment_representation(
        subfragment_representation, fragments, is_circular
    )
    return out_dseqrecord


def annotate_primer_binding_sites(
    input_dseqr: Dseqrecord, fragments: list[Dseqrecord]
) -> Dseqrecord:
    """Annotate the primer binding sites in a Dseqrecord."""
    fwd, _, rvs = fragments
    start_rvs = len(input_dseqr) - len(rvs)

    output_dseqr = copy.deepcopy(input_dseqr)
    output_dseqr.add_feature(
        x=0,
        y=len(fwd),
        type_="primer_bind",
        strand=1,
        label=[fwd.name],
        note=["sequence: " + str(fwd.seq)],
    )
    output_dseqr.add_feature(
        x=start_rvs,
        y=len(output_dseqr),
        type_="primer_bind",
        strand=-1,
        label=[rvs.name],
        note=["sequence: " + str(rvs.seq)],
    )
    return output_dseqr


def edge_representation2subfragment_representation(
    assembly: EdgeRepresentationAssembly, is_circular: bool
) -> SubFragmentRepresentationAssembly:
    """
    Turn this kind of edge representation fragment 1, fragment 2, right edge on 1, left edge on 2
    a = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b', 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    Into this: fragment 1, left edge on 1, right edge on 1
    b = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    """

    if is_circular:
        temp = list(assembly[-1:]) + list(assembly)
    else:
        temp = (
            [(None, assembly[0][0], None, None)]
            + list(assembly)
            + [(assembly[-1][1], None, None, None)]
        )
    edge_pairs = zip(temp, temp[1:])
    subfragment_representation = list()
    for (_u1, v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
        subfragment_representation.append((v1, start_location, end_location))

    return tuple(subfragment_representation)


def subfragment_representation2edge_representation(
    assembly: SubFragmentRepresentationAssembly, is_circular: bool
) -> EdgeRepresentationAssembly:
    """
    Turn this kind of subfragment representation fragment 1, left edge on 1, right edge on 1
    a = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    Into this: fragment 1, fragment 2, right edge on 1, left edge on 2
    b = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b' 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    """

    edge_representation = []

    # Iterate through the assembly pairwise to create the edge representation
    for i in range(len(assembly) - 1):
        frag1, left1, right1 = assembly[i]
        frag2, left2, right2 = assembly[i + 1]
        # Create the edge between the current and next fragment
        edge_representation.append((frag1, frag2, right1, left2))

    if is_circular:
        # Add the edge from the last fragment back to the first
        frag_last, left_last, right_last = assembly[-1]
        frag_first, left_first, right_first = assembly[0]
        edge_representation.append((frag_last, frag_first, right_last, left_first))

    return tuple(edge_representation)


def get_assembly_subfragments(
    fragments: list[Dseqrecord],
    subfragment_representation: SubFragmentRepresentationAssembly,
) -> list[Dseqrecord]:
    """From the fragment representation returned by edge_representation2subfragment_representation, get the subfragments that are joined together.

    Subfragments are the slices of the fragments that are joined together

    For example::

          --A--
        TACGTAAT
          --B--
         TCGTAACGA

        Gives: TACGTAA / CGTAACGA

    To reproduce::

        a = Dseqrecord('TACGTAAT')
        b = Dseqrecord('TCGTAACGA')
        f = Assembly([a, b], limit=5)
        a0 = f.get_linear_assemblies()[0]
        print(assembly2str(a0))
        a0_subfragment_rep =edge_representation2subfragment_representation(a0, False)
        for f in get_assembly_subfragments([a, b], a0_subfragment_rep):
            print(f.seq)

        # prints TACGTAA and CGTAACGA

    Subfragments: ``cccccgtatcgtgt``, ``atcgtgtactgtcatattc``
    """
    subfragments = list()
    for node, start_location, end_location in subfragment_representation:
        seq = (
            fragments[node - 1]
            if node > 0
            else fragments[-node - 1].reverse_complement()
        )
        subfragments.append(extract_subfragment(seq, start_location, end_location))
    return subfragments


def extract_subfragment(
    seq: Dseqrecord, start_location: Location | None, end_location: Location | None
) -> Dseqrecord:
    """Extract a subfragment from a sequence for an assembly, given the start and end locations of the subfragment."""

    if seq.circular and (start_location is None or end_location is None):
        raise ValueError(
            "Start and end locations cannot be None for circular sequences"
        )
        # This could be used to have consistent behaviour for circular sequences, where the start is arbitrary. However,
        # they should never get None, so this is not used.
        # if start_location is None:
        #     start_location = end_location
        # elif end_location is None:
        #     end_location = start_location

    start = 0 if start_location is None else location_boundaries(start_location)[0]
    end = None if end_location is None else location_boundaries(end_location)[1]

    # Special case, some of it could be handled by better Dseqrecord slicing in the future
    if seq.circular and locations_overlap(start_location, end_location, len(seq)):
        # The overhang is different for origin-spanning features, for instance
        # for a feature join{[12:13], [0:3]} in a sequence of length 13, the overhang
        # is -4, not 9
        ovhg = start - end if end > start else start - end - len(seq)
        # edge case
        if abs(ovhg) == len(seq):
            ovhg = 0
        dummy_cut = ((start, ovhg), None)
        open_seq = seq.apply_cut(dummy_cut, dummy_cut)
        return Dseqrecord(open_seq.seq.cast_to_ds(), features=open_seq.features)

    return seq[start:end]


def is_sublist(sublist: list, my_list: list, my_list_is_cyclic: bool = False) -> bool:
    """Returns True if argument sublist is a sublist of argument my_list (can be treated as cyclic), False otherwise.

    Examples
    --------
    >>> is_sublist([1, 2], [1, 2, 3], False)
    True
    >>> is_sublist([1, 2], [1, 3, 2], False)
    False

    # See the case here for cyclic lists
    >>> is_sublist([3, 1], [1, 2, 3], False)
    False
    >>> is_sublist([3, 1], [1, 2, 3], True)
    True
    """
    n = len(sublist)
    if my_list_is_cyclic:
        my_list = my_list + my_list
    for i in range(len(my_list) - n + 1):
        # Just in case tuples were passed
        if list(my_list[i : i + n]) == list(sublist):
            return True
    return False


def circular_permutation_min_abs(lst: list) -> list:
    """Returns the circular permutation of lst with the smallest absolute value first.

    Examples
    --------
    >>> circular_permutation_min_abs([1, 2, 3])
    [1, 2, 3]
    >>> circular_permutation_min_abs([3, 1, 2])
    [1, 2, 3]
    """
    min_abs_index = min(range(len(lst)), key=lambda i: abs(lst[i]))
    return lst[min_abs_index:] + lst[:min_abs_index]


class Assembly:
    """Assembly of a list of DNA fragments into linear or circular constructs.
    Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object. Several methods are available for analysis
    of overlapping sequences, graph construction and assembly.

    The assembly contains a directed graph, where nodes represent fragments and
    edges represent overlaps between fragments. :

    - The node keys are integers, representing the index of the fragment in the
      input list of fragments. The sign of the node key represents the orientation
      of the fragment, positive for forward orientation, negative for reverse orientation.
    - The edges contain the locations of the overlaps in the fragments. For an edge (u, v, key):
        - u and v are the nodes connected by the edge.
        - key is a string that represents the location of the overlap. In the format:
          'u[start:end](strand):v[start:end](strand)'.
        - Edges have a 'locations' attribute, which is a list of two FeatureLocation objects,
          representing the location of the overlap in the u and v fragment, respectively.
        - You can think of an edge as a representation of the join of two fragments.

    If fragment 1 and 2 share a subsequence of 6bp, [8:14] in fragment 1 and [1:7] in fragment 2,
    there will be 4 edges representing that overlap in the graph, for all possible
    orientations of the fragments (see add_edges_from_match for details):

    - ``(1, 2, '1[8:14]:2[1:7]')``
    - ``(2, 1, '2[1:7]:1[8:14]')``
    - ``(-1, -2, '-1[0:6]:-2[10:16]')``
    - ``(-2, -1, '-2[10:16]:-1[0:6]')``

    An assembly can be thought of as a tuple of graph edges, but instead of representing them with node indexes and keys, we represent them
    as u, v, locu, locv, where u and v are the nodes connected by the edge, and locu and locv are the locations of the overlap in the first
    and second fragment. Assemblies are then represented as:

    - Linear: ((1, 2, [8:14], [1:7]), (2, 3, [10:17], [1:8]))
    - Circular: ((1, 2, [8:14], [1:7]), (2, 3, [10:17], [1:8]), (3, 1, [12:17], [1:6]))

    Note that the first and last fragment are the same in a circular assembly.

    The following constrains are applied to remove duplicate assemblies:

    - Circular assemblies: the first subfragment is not reversed, and has the smallest index in the input fragment list.
      use_fragment_order is ignored.
    - Linear assemblies:
        - Using uid (see add_edges_from_match) to identify unique edges.

    Parameters
    ----------
    frags : list
        A list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered, this is passed as the third argument to the ``algorithm`` function.
        For certain algorithms, this might be ignored.
    algorithm : function, optional
        The algorithm used to determine the shared sequences. It's a function that takes two Dseqrecord objects as inputs,
        and will get passed the third argument (limit), that may or may not be used. It must return a list of overlaps
        (see common_sub_strings for an example).
    use_fragment_order : bool, optional
        It's set to True by default to reproduce legacy pydna behaviour: only assemblies that start with the first fragment and end with the last are considered.
        You should set it to False.
    use_all_fragments : bool, optional
        Constrain the assembly to use all fragments.


    Examples
    --------

    from assembly2 import Assembly, assembly2str
    from pydna.dseqrecord import Dseqrecord

    example_fragments = (
        Dseqrecord('AacgatCAtgctcc', name='a'),
        Dseqrecord('TtgctccTAAattctgc', name='b'),
        Dseqrecord('CattctgcGAGGacgatG', name='c'),
    )

    asm = Assembly(example_fragments, limit=5, use_fragment_order=False)
    print('Linear ===============')
    for assembly in asm.get_linear_assemblies():
        print(' ', assembly2str(assembly))
    print('Circular =============')
    for assembly in asm.get_circular_assemblies():
        print(' ', assembly2str(assembly))

    # Prints
    Linear ===============
        ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]')
        ('2[10:17]:3[1:8]', '3[12:17]:1[1:6]')
        ('3[12:17]:1[1:6]', '1[8:14]:2[1:7]')
        ('1[1:6]:3[12:17]',)
        ('2[1:7]:1[8:14]',)
        ('3[1:8]:2[10:17]',)
    Circular =============
        ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]', '3[12:17]:1[1:6]')

    """

    def __init__(
        self,
        frags: list[Dseqrecord],
        limit: int = 25,
        algorithm: AssemblyAlgorithmType = common_sub_strings,
        use_fragment_order: bool = True,
        use_all_fragments: bool = False,
    ):

        # TODO: allow for the same fragment to be included more than once?
        self.G = nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_nodes_from((i + 1, {"seq": f}) for (i, f) in enumerate(frags))
        self.G.add_nodes_from(
            (-(i + 1), {"seq": f.reverse_complement()}) for (i, f) in enumerate(frags)
        )

        # Iterate over all possible combinations of fragments
        fragment_pairs = itertools.combinations(
            filter(lambda x: x > 0, self.G.nodes), 2
        )
        for i, j in fragment_pairs:
            # All the relative orientations of the fragments in the pair
            for u, v in itertools.product([i, -i], [j, -j]):
                u_seq = self.G.nodes[u]["seq"]
                v_seq = self.G.nodes[v]["seq"]
                matches = algorithm(u_seq, v_seq, limit)
                for match in matches:
                    self.add_edges_from_match(match, u, v, u_seq, v_seq)

        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order
        self.use_all_fragments = use_all_fragments

        return

    @classmethod
    def assembly_is_valid(
        cls,
        fragments: list[Dseqrecord | Primer],
        assembly: EdgeRepresentationAssembly,
        is_circular: bool,
        use_all_fragments: bool,
        is_insertion: bool = False,
    ) -> bool:
        """
        Returns True if the assembly is valid, False otherwise. See function comments for conditions tested.
        """
        if is_circular is None:
            return False

        # Linear assemblies may get begin-1-end, begin-2-end, these are removed here.
        if len(assembly) == 0:
            return False

        # Topology check -> Circular sequences cannot be first or last in a linear assembly.
        # For example, let's imagine aACGTc (linear) and gACGTc (circular).
        # It should not be possible to join them into a linear assembly. It's similar if we
        # think of a restriction-ligation assembly, example: aGAATTCc (linear) and gGAATTCc
        # (circular).
        # A linear product can be generated where the circular molecule is cut open, and one end
        # it joins the linear molecule and on the other it's free, but for now it's not a
        # relevant product and it's excluded.
        first_fragment = fragments[abs(assembly[0][0]) - 1]
        last_fragment = fragments[abs(assembly[-1][1]) - 1]
        if not is_circular and (
            isinstance(first_fragment, Dseqrecord)
            and first_fragment.circular
            or (isinstance(last_fragment, Dseqrecord) and last_fragment.circular)
        ):
            return False

        if use_all_fragments and len(fragments) != len(
            set(flatten(map(abs, e[:2]) for e in assembly))
        ):
            return False

        # Here we check whether subsequent pairs of fragments are compatible, for instance:
        # Compatible (overlap of 1 and 2 occurs before overlap of 2 and 3):
        # (1,2,[2:9],[0:7]), (2,3,[12:19],[0:7])
        #    -- A --
        # 1 gtatcgtgt     -- B --
        # 2   atcgtgtactgtcatattc
        # 3               catattcaa
        # Incompatible (overlap of 1 and 2 occurs after overlap of 2 and 3):
        # (1,2,[2:9],[13:20]), (2,3,[0:7],[0:7])
        #                 -- A --
        #  1 -- B --    gtatcgtgt
        #  2 catattcccccccatcgtgtactgt
        #  3 catattcaa
        # Redundant: overlap of 1 and 2 ends at the same spot as overlap of 2 and 3
        # (1,2,[2:9],[1:8]), (2,3,[0:8],[0:8])
        #    -- A --
        #  gtatcgtgt
        #   catcgtgtactgtcatattc
        #   catcgtgtactgtcatattc
        #   -- B ---
        if is_circular:
            # In a circular assembly, first and last fragment must be the same
            if assembly[0][0] != assembly[-1][1]:
                return False
            edge_pairs = zip(assembly, assembly[1:] + assembly[:1])
        else:
            edge_pairs = zip(assembly, assembly[1:])

        for (_u1, v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
            # Incompatible as described in figure above
            fragment = fragments[abs(v1) - 1]
            if (
                isinstance(fragment, Primer) or not fragment.circular
            ) and location_boundaries(start_location)[1] >= location_boundaries(
                end_location
            )[
                1
            ]:
                return False

        # Fragments are used only once
        nodes_used = [
            f[0]
            for f in edge_representation2subfragment_representation(
                assembly, is_circular or is_insertion
            )
        ]
        if len(nodes_used) != len(set(map(abs, nodes_used))):
            return False

        return True

    def add_edges_from_match(
        self,
        match: SequenceOverlap,
        u: int,
        v: int,
        first: Dseqrecord,
        secnd: Dseqrecord,
    ):
        """Add edges to the graph from a match returned by the ``algorithm`` function (see pydna.common_substrings). For
        format of edges (see documentation of the Assembly class).

        Matches are directional, because not all ``algorithm`` functions return the same match for (u,v) and (v,u). For example,
        homologous recombination does but sticky end ligation does not. The function returns two edges:

        - Fragments in the orientation they were passed, with locations of the match (u, v, loc_u, loc_v)
        - Reverse complement of the fragments with inverted order, with flipped locations (-v, -u, flip(loc_v), flip(loc_u))/

        """
        x_start, y_start, length = match
        if length == 0:
            # Edge case, blunt ligation
            locs = [SimpleLocation(x_start, x_start), SimpleLocation(y_start, y_start)]
        else:
            # We use shift_location with 0 to wrap origin-spanning features
            locs = [
                shift_location(
                    SimpleLocation(x_start, x_start + length), 0, len(first)
                ),
                shift_location(
                    SimpleLocation(y_start, y_start + length), 0, len(secnd)
                ),
            ]

        # Flip the locations to get the reverse complement
        rc_locs = [locs[0]._flip(len(first)), locs[1]._flip(len(secnd))]

        # Unique id that identifies the edge in either orientation
        uid = f"{u}{locs[0]}:{v}{locs[1]}"

        combinations = (
            (u, v, locs),
            (-v, -u, rc_locs[::-1]),
        )

        for u, v, l in combinations:
            self.G.add_edge(u, v, f"{u}{l[0]}:{v}{l[1]}", locations=l, uid=uid)

    def format_assembly_edge(
        self, graph_edge: tuple[int, int, str]
    ) -> AssemblyEdgeType:
        """Go from the (u, v, key) to the (u, v, locu, locv) format."""
        u, v, key = graph_edge
        locu, locv = self.G.get_edge_data(u, v, key)["locations"]
        return u, v, locu, locv

    def get_linear_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        """Get linear assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid). Subassemblies are removed (see remove_subassemblies).
        """

        # Copy the graph since we will add the begin and end mock nodes
        G = nx.MultiDiGraph(self.G)
        G.add_nodes_from(["begin", "end"])

        if self.use_fragment_order:
            # Path must start with the first fragment and end with the last
            G.add_edge("begin", 1)
            G.add_edge("begin", -1)
            G.add_edge(len(self.fragments), "end")
            G.add_edge(-len(self.fragments), "end")
        else:
            for node in filter(lambda x: type(x) is int, G.nodes):
                G.add_edge("begin", node)
                G.add_edge(node, "end")

        unique_linear_paths = self.get_unique_linear_paths(G)
        possible_assemblies = self.get_possible_assembly_number(unique_linear_paths)
        if possible_assemblies > max_assemblies:
            raise ValueError(
                f"Too many assemblies ({possible_assemblies} pre-validation) to assemble"
            )

        assemblies = sum(
            map(lambda x: self.node_path2assembly_list(x, False), unique_linear_paths),
            [],
        )

        out = [
            a
            for a in assemblies
            if self.assembly_is_valid(self.fragments, a, False, self.use_all_fragments)
        ]
        if only_adjacent_edges:
            out = [a for a in out if self.assembly_uses_only_adjacent_edges(a, False)]
        return remove_subassemblies(out)

    def node_path2assembly_list(
        self, cycle: list[int], circular: bool
    ) -> list[EdgeRepresentationAssembly]:
        """Convert a node path in the format [1, 2, 3] (as returned by networkx.cycles.simple_cycles) to a list of all
          possible assemblies.

        There may be multiple assemblies for a given node path, if there are several edges connecting two nodes,
        for example two overlaps between 1 and 2, and single overlap between 2 and 3 should return 3 assemblies.
        """
        combine = list()
        pairing = (
            zip(cycle, cycle[1:] + cycle[:1]) if circular else zip(cycle, cycle[1:])
        )
        for u, v in pairing:
            combine.append([(u, v, key) for key in self.G[u][v]])
        return [
            tuple(map(self.format_assembly_edge, x))
            for x in itertools.product(*combine)
        ]

    def get_unique_linear_paths(
        self, G_with_begin_end: nx.MultiDiGraph, max_paths=10000
    ) -> list[list[int]]:
        """Get unique linear paths from the graph, removing those that contain the same node twice."""
        # We remove the begin and end nodes, and get all paths without edges
        # e.g. we will get [1, 2, 3] only once, even if multiple edges connect
        # 1 and 2 or 2 and 3, by converting to DiGraph.

        # Cutoff has a different meaning of what one would expect, see https://github.com/networkx/networkx/issues/2762
        node_paths = [
            x[1:-1]
            for x in limit_iterator(
                nx.all_simple_paths(
                    nx.DiGraph(G_with_begin_end),
                    "begin",
                    "end",
                    cutoff=(len(self.fragments) + 1),
                ),
                max_paths,
            )
        ]

        # Remove those that contain the same node twice
        node_paths = [x for x in node_paths if len(x) == len(set(map(abs, x)))]

        if self.use_all_fragments:
            node_paths = [x for x in node_paths if len(x) == len(self.fragments)]

        # For each path, we check if there are reverse complement duplicates
        # See: https://github.com/manulera/OpenCloning_backend/issues/160
        unique_node_paths = list()
        for p in node_paths:
            if [-x for x in p[::-1]] not in unique_node_paths:
                unique_node_paths.append(p)

        return unique_node_paths

    def get_possible_assembly_number(self, paths: list[list[int]]) -> int:
        """
        Get the number of possible assemblies from a list of node paths. Basically, for each path
        passed as a list of integers / nodes, we calculate the number of paths possible connecting
        the nodes in that order, given the graph (all the edges connecting them).
        """
        possibilities = 0
        for path in paths:
            this_path = 1
            for u, v in zip(path, path[1:]):
                if v in self.G[u]:
                    this_path *= len(self.G[u][v])
            possibilities += this_path
        return possibilities

    def get_circular_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        """Get circular assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid)."""
        # The constrain of circular sequence is that the first node is the fragment with the smallest index in its initial orientation,
        # this is ensured by the circular_permutation_min_abs function + the filter below
        sorted_cycles = map(
            circular_permutation_min_abs,
            limit_iterator(
                nx.cycles.simple_cycles(self.G, length_bound=len(self.fragments)),
                10000,
            ),
        )
        sorted_cycles = filter(lambda x: x[0] > 0, sorted_cycles)
        # cycles.simple_cycles returns lists [1,2,3] not assemblies, see self.cycle2circular_assemblies

        # We apply constrains already here because sometimes the combinatorial explosion is too large
        if self.use_all_fragments:
            sorted_cycles = [c for c in sorted_cycles if len(c) == len(self.fragments)]

        # Remove cycles with duplicates
        sorted_cycles = [c for c in sorted_cycles if len(c) == len(set(map(abs, c)))]
        possible_assembly_number = self.get_possible_assembly_number(
            [c + c[:1] for c in sorted_cycles]
        )
        if possible_assembly_number > max_assemblies:
            raise ValueError(
                f"Too many assemblies ({possible_assembly_number} pre-validation) to assemble"
            )

        assemblies = sum(
            map(lambda x: self.node_path2assembly_list(x, True), sorted_cycles), []
        )

        out = [
            a
            for a in assemblies
            if self.assembly_is_valid(self.fragments, a, True, self.use_all_fragments)
        ]
        if only_adjacent_edges:
            out = [a for a in out if self.assembly_uses_only_adjacent_edges(a, True)]
        return out

    def format_insertion_assembly(
        self, assembly: EdgeRepresentationAssembly
    ) -> EdgeRepresentationAssembly | None:
        """Sorts the fragment representing a cycle so that they represent an insertion assembly if possible,
        else returns None.

        Here we check if one of the joins between fragments represents the edges of an insertion assembly
        The fragment must be linear, and the join must be as indicated below

        ::

            --------         -------           Fragment 1
                ||            ||
                xxxxxxxx      ||               Fragment 2
                      ||      ||
                      oooooooooo               Fragment 3

        The above example will be [(1, 2, [4:6], [0:2]), (2, 3, [6:8], [0:2]), (3, 1, [8:10], [9:11)])]

        These could be returned in any order by simple_cycles, so we sort the edges so that the first
        and last ``u`` and ``v`` match the fragment that gets the insertion (1 in the example above).
        """
        edge_pair_index = list()

        # Pair edges with one another
        for i, ((_u1, v1, _, end_location), (_u2, _v2, start_location, _)) in enumerate(
            zip(assembly, assembly[1:] + assembly[:1])
        ):
            fragment = self.fragments[abs(v1) - 1]
            # Find the pair of edges that should be last and first  ((3, 1, [8:10], [9:11)]), (1, 2, [4:6], [0:2]) in
            # the example above. Only one of the pairs of edges should satisfy this condition for the topology to make sense.
            left_of_insertion = location_boundaries(start_location)[0]
            right_of_insertion = location_boundaries(end_location)[0]
            if not fragment.circular and (
                right_of_insertion >= left_of_insertion
                # The below condition is for single-site integration.
                # The reason to use locations_overlap instead of equality is because the location might extend
                # left of right. For example, let's take ACCGGTTT as homology arm for an integration:
                #
                # insert aaACCGGTTTccACCGGTTTtt
                # genome aaACCGGTTTtt
                #
                # The locations of homology on the genome are [0:10] and [2:12], so not identical
                # but they overlap.
                or locations_overlap(start_location, end_location, len(fragment))
            ):
                edge_pair_index.append(i)

        if len(edge_pair_index) != 1:
            return None

        shift_by = (edge_pair_index[0] + 1) % len(assembly)
        return assembly[shift_by:] + assembly[:shift_by]

    def format_insertion_assembly_edge_case(
        self, assembly: EdgeRepresentationAssembly
    ) -> EdgeRepresentationAssembly:
        """
        Edge case from https://github.com/manulera/OpenCloning_backend/issues/329
        """
        same_assembly = assembly[:]

        if len(assembly) != 2:
            return same_assembly
        ((f1, f2, loc_f1_1, loc_f2_1), (_f2, _f1, loc_f2_2, loc_f1_2)) = assembly

        if f1 != _f1 or _f2 != f2:
            return same_assembly

        if loc_f2_1 == loc_f2_2 or loc_f1_2 == loc_f1_1:
            return same_assembly

        fragment1 = self.fragments[abs(f1) - 1]
        fragment2 = self.fragments[abs(f2) - 1]

        if not locations_overlap(
            loc_f1_1, loc_f1_2, len(fragment1)
        ) or not locations_overlap(loc_f2_2, loc_f2_1, len(fragment2)):
            return same_assembly

        # Sort to make compatible with insertion assembly
        if location_boundaries(loc_f1_1)[0] > location_boundaries(loc_f1_2)[0]:
            new_assembly = same_assembly[::-1]
        else:
            new_assembly = same_assembly[:]

        ((f1, f2, loc_f1_1, loc_f2_1), (_f2, _f1, loc_f2_2, loc_f1_2)) = new_assembly

        fragment1 = self.fragments[abs(f1) - 1]
        if fragment1.circular:
            return same_assembly
        fragment2 = self.fragments[abs(f2) - 1]

        # Extract boundaries
        f2_1_start, _ = location_boundaries(loc_f2_1)
        f2_2_start, f2_2_end = location_boundaries(loc_f2_2)
        f1_1_start, _ = location_boundaries(loc_f1_1)
        f1_2_start, f1_2_end = location_boundaries(loc_f1_2)

        overlap_diff = len(fragment1[f1_1_start:f1_2_end]) - len(
            fragment2[f2_1_start:f2_2_end]
        )

        # Safeguard
        if overlap_diff == 0:  # pragma: no cover
            raise AssertionError("Overlap is 0")

        if overlap_diff > 0:
            new_loc_f1_1 = create_location(
                f1_1_start, f1_2_start - overlap_diff, len(fragment1)
            )
            new_loc_f2_1 = create_location(f2_1_start, f2_2_start, len(fragment2))
        else:
            new_loc_f2_1 = create_location(
                f2_1_start, f2_2_start + overlap_diff, len(fragment2)
            )
            new_loc_f1_1 = create_location(f1_1_start, f1_2_start, len(fragment1))

        new_assembly = [
            (f1, f2, new_loc_f1_1, new_loc_f2_1),
            new_assembly[1],
        ]

        return new_assembly

    def get_insertion_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        """Assemblies that represent the insertion of a fragment or series of fragment inside a linear construct. For instance,
        digesting CCCCGAATTCCCCGAATTC with EcoRI and inserting the fragment with two overhangs into the EcoRI site of AAAGAATTCAAA.
        This is not so much meant for the use-case of linear fragments that represent actual linear fragments, but for linear
        fragments that represent a genome region. This can then be used to simulate homologous recombination.
        """
        if only_adjacent_edges:
            raise NotImplementedError(
                "only_adjacent_edges not implemented for insertion assemblies"
            )

        cycles = limit_iterator(nx.cycles.simple_cycles(self.G), 10000)

        # We apply constrains already here because sometimes the combinatorial explosion is too large
        if self.use_all_fragments:
            cycles = [c for c in cycles if len(c) == len(self.fragments)]

        # Remove cycles with duplicates
        cycles = [c for c in cycles if len(c) == len(set(map(abs, c)))]

        possible_assembly_number = self.get_possible_assembly_number(
            [c + c[:1] for c in cycles]
        )

        if possible_assembly_number > max_assemblies:
            raise ValueError(
                f"Too many assemblies ({possible_assembly_number} pre-validation) to assemble"
            )

        # We find cycles first
        iterator = limit_iterator(nx.cycles.simple_cycles(self.G), 10000)
        assemblies = sum(
            map(lambda x: self.node_path2assembly_list(x, True), iterator), []
        )
        # We format the edge case
        assemblies = [self.format_insertion_assembly_edge_case(a) for a in assemblies]
        # We select those that contain exactly only one suitable edge
        assemblies = [
            b
            for a in assemblies
            if (b := self.format_insertion_assembly(a)) is not None
        ]
        # First fragment should be in the + orientation
        assemblies = list(filter(lambda x: x[0][0] > 0, assemblies))
        return [
            a
            for a in assemblies
            if self.assembly_is_valid(
                self.fragments, a, False, self.use_all_fragments, is_insertion=True
            )
        ]

    def assemble_linear(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[Dseqrecord]:
        """Assemble linear constructs, from assemblies returned by self.get_linear_assemblies."""
        assemblies = self.get_linear_assemblies(only_adjacent_edges, max_assemblies)
        return [assemble(self.fragments, a) for a in assemblies]

    def assemble_circular(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[Dseqrecord]:
        """Assemble circular constructs, from assemblies returned by self.get_circular_assemblies."""
        assemblies = self.get_circular_assemblies(only_adjacent_edges, max_assemblies)
        return [assemble(self.fragments, a) for a in assemblies]

    def assemble_insertion(self, only_adjacent_edges: bool = False) -> list[Dseqrecord]:
        """Assemble insertion constructs, from assemblies returned by self.get_insertion_assemblies."""
        assemblies = self.get_insertion_assemblies(only_adjacent_edges)
        return [assemble(self.fragments, a, is_insertion=True) for a in assemblies]

    def get_locations_on_fragments(self) -> dict[int, dict[str, list[Location]]]:
        """Get a dictionary where the keys are the nodes in the graph, and the values are dictionaries with keys
        ``left``, ``right``, containing (for each fragment) the locations where the fragment is joined to another fragment on its left
        and right side. The values in ``left`` and ``right`` are often the same, except in restriction-ligation with partial overlap enabled,
        where we can end up with a situation like this:

        GGTCTCCCCAATT and aGGTCTCCAACCAA as fragments

        # Partial overlap in assembly 1[9:11]:2[8:10]
        GGTCTCCxxAACCAA
        CCAGAGGGGTTxxTT

        # Partial overlap in 2[10:12]:1[7:9]
        aGGTCTCCxxCCAATT
        tCCAGAGGTTGGxxAA

        Would return::

            {
                1: {'left': [7:9], 'right': [9:11]},
                2: {'left': [8:10], 'right': [10:12]},
                -1: {'left': [2:4], 'right': [4:6]},
                -2: {'left': [2:4], 'right': [4:6]}
            }

        """

        locations_on_fragments = dict()
        for node in self.G.nodes:
            this_dict = {"left": list(), "right": list()}
            for edge in self.G.edges(data=True):
                for i, key in enumerate(["right", "left"]):
                    if edge[i] == node:
                        edge_location = edge[2]["locations"][i]
                        if edge_location not in this_dict[key]:
                            this_dict[key].append(edge_location)
            this_dict["left"] = sorted(
                this_dict["left"], key=lambda x: location_boundaries(x)[0]
            )
            this_dict["right"] = sorted(
                this_dict["right"], key=lambda x: location_boundaries(x)[0]
            )
            locations_on_fragments[node] = this_dict

        return locations_on_fragments

    def assembly_uses_only_adjacent_edges(self, assembly, is_circular: bool) -> bool:
        """
        Check whether only adjacent edges within each fragment are used in the assembly. This is useful to check if a cut and ligate assembly is valid,
        and prevent including partially digested fragments. For example, imagine the following fragment being an input for a digestion
        and ligation assembly, where the enzyme cuts at the sites indicated by the vertical lines:

        ::

                     x       y       z
              -------|-------|-------|---------

        We would only want assemblies that contain subfragments start-x, x-y, y-z, z-end, and not start-x, y-end, for instance.
        The latter would indicate that the fragment was partially digested.
        """

        locations_on_fragments = self.get_locations_on_fragments()
        for node in locations_on_fragments:
            fragment_len = len(self.fragments[abs(node) - 1])
            for side in ["left", "right"]:
                locations_on_fragments[node][side] = gather_overlapping_locations(
                    locations_on_fragments[node][side], fragment_len
                )

        allowed_location_pairs = dict()
        for node in locations_on_fragments:
            if not is_circular:
                # We add the existing ends of the fragment
                left = [(None,)] + locations_on_fragments[node]["left"]
                right = locations_on_fragments[node]["right"] + [(None,)]

            else:
                # For circular assemblies, we add the first location at the end
                # to allow for the last edge to be used
                left = locations_on_fragments[node]["left"]
                right = (
                    locations_on_fragments[node]["right"][1:]
                    + locations_on_fragments[node]["right"][:1]
                )

            pairs = list()
            for pair in zip(left, right):
                pairs += list(itertools.product(*pair))
            allowed_location_pairs[node] = pairs

        fragment_assembly = edge_representation2subfragment_representation(
            assembly, is_circular
        )
        for node, start_location, end_location in fragment_assembly:
            if (start_location, end_location) not in allowed_location_pairs[node]:
                return False
        return True

    def __repr__(self):
        # https://pyformat.info
        return ps(
            "Assembly\n"
            "fragments..: {sequences}\n"
            "limit(bp)..: {limit}\n"
            "G.nodes....: {nodes}\n"
            "algorithm..: {al}".format(
                sequences=" ".join("{}bp".format(len(x)) for x in self.fragments),
                limit=self.limit,
                nodes=self.G.order(),
                al=self.algorithm.__name__,
            )
        )


class PCRAssembly(Assembly):
    """
    An assembly that represents a PCR, where ``fragments`` is a list of primer, template, primer (in that order).
    It always uses the ``primer_template_overlap`` algorithm and accepts the ``mismatches`` argument to indicate
    the number of mismatches allowed in the overlap. Only supports substitution mismatches, not indels.
    """

    def __init__(self, frags: list[Dseqrecord | Primer], limit=25, mismatches=0):

        value_error = ValueError(
            "PCRAssembly assembly must be initialised with a list/tuple of primer, template, primer"
        )
        if len(frags) != 3:
            raise value_error

        # Validate the inputs: should be a series of primer, template, primer
        wrong_fragment_class = (
            not isinstance(frags[0], Primer),
            isinstance(frags[1], Primer),
            not isinstance(frags[2], Primer),
        )
        if any(wrong_fragment_class):
            raise value_error

        # TODO: allow for the same fragment to be included more than once?
        self.G = nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_nodes_from((i + 1, {"seq": f}) for (i, f) in enumerate(frags))
        self.G.add_nodes_from(
            (-(i + 1), {"seq": f.reverse_complement()}) for (i, f) in enumerate(frags)
        )

        pairs = list()
        primer_ids = list()
        for i in range(0, len(frags), 3):
            # primer, template, primer
            p1, t, p2 = (i + 1, i + 2, i + 3)
            primer_ids += [p1, p2]
            pairs += list(itertools.product([p1, p2], [t, -t]))
            pairs += list(itertools.product([t, -t], [-p1, -p2]))

        for u, v in pairs:
            u_seq = self.G.nodes[u]["seq"]
            v_seq = self.G.nodes[v]["seq"]
            matches = primer_template_overlap(u_seq, v_seq, limit, mismatches)
            for match in matches:
                self.add_edges_from_match(match, u, v, u_seq, v_seq)

        # These two are constrained
        self.use_fragment_order = False
        self.use_all_fragments = True

        self.fragments = frags
        self.limit = limit
        self.algorithm = primer_template_overlap

        return

    def get_linear_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        if only_adjacent_edges:
            raise NotImplementedError(
                "only_adjacent_edges not implemented for PCR assemblies"
            )

        return super().get_linear_assemblies(max_assemblies=max_assemblies)

    def get_circular_assemblies(self, only_adjacent_edges: bool = False):
        raise NotImplementedError(
            "get_circular_assemblies not implemented for PCR assemblies"
        )

    def get_insertion_assemblies(self, only_adjacent_edges: bool = False):
        raise NotImplementedError(
            "get_insertion_assemblies not implemented for PCR assemblies"
        )

    def assemble_linear(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[Dseqrecord]:
        """
        Overrides the parent method to ensure that the 5' of the crick strand of the product matches the
        sequence of the reverse primer. This is important when using primers with dUTP (for USER cloning).
        """
        results = super().assemble_linear(only_adjacent_edges, max_assemblies)
        for result in results:
            rp = self.fragments[2]
            result.seq = result.seq[: -len(rp)] + Dseq(str(rp.seq.rc()))
        return results


class SingleFragmentAssembly(Assembly):
    """
    An assembly that represents the circularisation or splicing of a single fragment.
    """

    def __init__(self, frags: [Dseqrecord], limit=25, algorithm=common_sub_strings):

        if len(frags) != 1:
            raise ValueError(
                "SingleFragmentAssembly assembly must be initialised with a single fragment"
            )
        # TODO: allow for the same fragment to be included more than once?
        self.G = nx.MultiDiGraph()
        frag = frags[0]
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_node(1, seq=frag)

        matches = algorithm(frag, frag, limit)
        for match in matches:
            self.add_edges_from_match(match, 1, 1, frag, frag)

        # To avoid duplicated outputs
        self.G.remove_edges_from([(-1, -1)])

        # These two are constrained
        self.use_fragment_order = True
        self.use_all_fragments = True

        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm

        return

    def get_circular_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        # We don't want the same location twice
        assemblies = filter(
            lambda x: x[0][2] != x[0][3],
            super().get_circular_assemblies(only_adjacent_edges, max_assemblies),
        )
        return [
            a
            for a in assemblies
            if self.assembly_is_valid(self.fragments, a, True, self.use_all_fragments)
        ]

    def get_insertion_assemblies(
        self, only_adjacent_edges: bool = False, max_assemblies: int = 50
    ) -> list[EdgeRepresentationAssembly]:
        """This could be renamed splicing assembly, but the essence is similar"""

        if only_adjacent_edges:
            raise NotImplementedError(
                "only_adjacent_edges not implemented for insertion assemblies"
            )

        def splicing_assembly_filter(x):
            # We don't want the same location twice
            if x[0][2] == x[0][3]:
                return False
            # We don't want to get overlap only (e.g. GAATTCcatGAATTC giving GAATTC)
            left_start, _ = location_boundaries(x[0][2])
            _, right_end = location_boundaries(x[0][3])
            if left_start == 0 and right_end == len(self.fragments[0]):
                return False
            return True

        # We don't want the same location twice
        assemblies = filter(
            splicing_assembly_filter,
            super().get_insertion_assemblies(max_assemblies=max_assemblies),
        )
        return [
            a
            for a in assemblies
            if self.assembly_is_valid(
                self.fragments, a, False, self.use_all_fragments, is_insertion=True
            )
        ]

    def get_linear_assemblies(self):
        raise NotImplementedError("Linear assembly does not make sense")


def common_function_assembly_products(
    frags: list[Dseqrecord],
    limit: int | None,
    algorithm: Callable,
    circular_only: bool,
    filter_results_function: Callable | None = None,
    only_adjacent_edges: bool = False,
) -> list[Dseqrecord]:
    """Common function to avoid code duplication. Could be simplified further
    once SingleFragmentAssembly and Assembly are merged.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    limit : int or None
        Minimum overlap length required, or None if not applicable
    algorithm : Callable
        Function that determines valid overlaps between fragments
    circular_only : bool
        If True, only return circular assemblies
    filter_results_function : Callable or None
        Function that filters the results
    only_adjacent_edges : bool
        If True, only return assemblies that use only adjacent edges

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """
    if len(frags) == 1:
        asm = SingleFragmentAssembly(frags, limit, algorithm)
    else:
        asm = Assembly(
            frags, limit, algorithm, use_fragment_order=False, use_all_fragments=True
        )
    output_assemblies = asm.get_circular_assemblies(only_adjacent_edges)
    if not circular_only and len(frags) > 1:
        output_assemblies += filter_linear_subassemblies(
            asm.get_linear_assemblies(only_adjacent_edges), output_assemblies, frags
        )
    if not circular_only and len(frags) == 1:
        output_assemblies += asm.get_insertion_assemblies()

    if filter_results_function:
        output_assemblies = [a for a in output_assemblies if filter_results_function(a)]

    return [assemble(frags, a) for a in output_assemblies]


def _recast_sources(
    products: list[Dseqrecord], source_cls, **extra_fields
) -> list[Dseqrecord]:
    """Recast the `source` of each product to `source_cls` with optional extras.

    This avoids repeating the same for-loop across many assembly functions.
    """
    for prod in products:
        prod.source = source_cls(
            **prod.source.to_unserialized_dict(),
            **extra_fields,
        )
    return products


def gibson_assembly(
    frags: list[Dseqrecord], limit: int = 25, circular_only: bool = False
) -> list[Dseqrecord]:
    """Returns the products for Gibson assembly.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    limit : int, optional
        Minimum overlap length required, by default 25
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """

    products = common_function_assembly_products(
        frags, limit, gibson_overlap, circular_only
    )
    return _recast_sources(products, GibsonAssemblySource)


def in_fusion_assembly(
    frags: list[Dseqrecord], limit: int = 25, circular_only: bool = False
) -> list[Dseqrecord]:
    """Returns the products for in-fusion assembly. This is the same as Gibson
    assembly, but with a different name.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    limit : int, optional
        Minimum overlap length required, by default 25
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """

    products = common_function_assembly_products(
        frags, limit, in_fusion_overlap, circular_only
    )
    return _recast_sources(products, InFusionSource)


def fusion_pcr_assembly(
    frags: list[Dseqrecord], limit: int = 25, circular_only: bool = False
) -> list[Dseqrecord]:
    """Returns the products for fusion PCR assembly. This is the same as Gibson
    assembly, but with a different name.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    limit : int, optional
        Minimum overlap length required, by default 25
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """
    products = common_function_assembly_products(
        frags, limit, pcr_fusion_overlap, circular_only
    )
    return _recast_sources(products, OverlapExtensionPCRLigationSource)


def in_vivo_assembly(
    frags: list[Dseqrecord], limit: int = 25, circular_only: bool = False
) -> list[Dseqrecord]:
    """Returns the products for in vivo assembly (IVA), which relies on homologous recombination between the fragments.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    limit : int, optional
        Minimum overlap length required, by default 25
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """
    products = common_function_assembly_products(
        frags, limit, common_sub_strings, circular_only
    )
    return _recast_sources(products, InVivoAssemblySource)


def restriction_ligation_assembly(
    frags: list[Dseqrecord],
    enzymes: list["AbstractCut"],
    allow_blunt: bool = True,
    circular_only: bool = False,
) -> list[Dseqrecord]:
    """Returns the products for restriction ligation assembly:

    - Finds cutsites in the fragments
    - Finds all products that could be assembled by ligating the fragments based on those cutsites
    - Will NOT return products that combine an existing end with an end generated by the same enzyme (see example below)

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    enzymes : list[AbstractCut]
        List of restriction enzymes to use
    allow_blunt : bool, optional
        If True, allow blunt end ligations, by default True
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules

    Examples
    --------
    In the example below, we plan to assemble a plasmid from a backbone and an insert, using the EcoRI and SalI enzymes.
    Note how 2 circular products are returned, one contains the insert (``acgt``)
    and the desired part of the backbone (``cccccc``), the other contains the
    reversed insert (``tgga``) and the cut-out part of the backbone (``aaa``).

    >>> from pydna.assembly2 import restriction_ligation_assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from Bio.Restriction import EcoRI, SalI
    >>> backbone = Dseqrecord("cccGAATTCaaaGTCGACccc", circular=True)
    >>> insert = Dseqrecord("ggGAATTCaggtGTCGACgg")
    >>> products = restriction_ligation_assembly([backbone, insert], [EcoRI, SalI], circular_only=True)
    >>> products[0].seq
    Dseq(o22)
    TCGACccccccGAATTCaggtG
    AGCTGggggggCTTAAGtccaC
    >>> products[1].seq
    Dseq(o19)
    AATTCaaaGTCGACacctG
    TTAAGtttCAGCTGtggaC

    Note that passing a pre-cut fragment will not work.

    >>> restriction_products = insert.cut([EcoRI, SalI])
    >>> cut_insert = restriction_products[1]
    >>> restriction_ligation_assembly([backbone, cut_insert], [EcoRI, SalI], circular_only=True)
    []

    It also works with a single fragment, for circularization:

    >>> seq = Dseqrecord("GAATTCaaaGAATTC")
    >>> products =restriction_ligation_assembly([seq], [EcoRI])
    >>> products[0].seq
    Dseq(o9)
    AATTCaaaG
    TTAAGtttC
    """

    def algorithm_fn(x, y, _l):
        # By default, we allow blunt ends
        return restriction_ligation_overlap(x, y, enzymes, False, allow_blunt)

    products = common_function_assembly_products(
        frags, None, algorithm_fn, circular_only, only_adjacent_edges=True
    )
    return _recast_sources(
        products, RestrictionAndLigationSource, restriction_enzymes=enzymes
    )


def golden_gate_assembly(
    frags: list[Dseqrecord],
    enzymes: list["AbstractCut"],
    allow_blunt: bool = True,
    circular_only: bool = False,
) -> list[Dseqrecord]:
    """Returns the products for Golden Gate assembly. This is the same as
    restriction ligation assembly, but with a different name. Check the documentation
    for ``restriction_ligation_assembly`` for more details.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    enzymes : list[AbstractCut]
        List of restriction enzymes to use
    allow_blunt : bool, optional
        If True, allow blunt end ligations, by default True
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules

    Examples
    --------
    See the example for ``restriction_ligation_assembly``.
    """
    return restriction_ligation_assembly(frags, enzymes, allow_blunt, circular_only)


def ligation_assembly(
    frags: list[Dseqrecord],
    allow_blunt: bool = False,
    allow_partial_overlap: bool = False,
    circular_only: bool = False,
) -> list[Dseqrecord]:
    """Returns the products for ligation assembly, as inputs pass the fragments (digested if needed) that
    will be ligated.

    For most cases, you probably should use ``restriction_ligation_assembly`` instead.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    allow_blunt : bool, optional
        If True, allow blunt end ligations, by default False
    allow_partial_overlap : bool, optional
        If True, allow partial overlaps between sticky ends, by default False
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules


    Examples
    --------
    In the example below, we plan to assemble a plasmid from a backbone and an insert,
    using the EcoRI enzyme. The insert and insertion site in the backbone are flanked by
    EcoRI sites, so there are two possible products depending on the orientation of the insert.

    >>> from pydna.assembly2 import ligation_assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from Bio.Restriction import EcoRI
    >>> backbone = Dseqrecord("cccGAATTCaaaGAATTCccc", circular=True)
    >>> backbone_cut = backbone.cut(EcoRI)[1]
    >>> insert = Dseqrecord("ggGAATTCaggtGAATTCgg")
    >>> insert_cut = insert.cut(EcoRI)[1]
    >>> products = ligation_assembly([backbone_cut, insert_cut])
    >>> products[0].seq
    Dseq(o22)
    AATTCccccccGAATTCaggtG
    TTAAGggggggCTTAAGtccaC
    >>> products[1].seq
    Dseq(o22)
    AATTCccccccGAATTCacctG
    TTAAGggggggCTTAAGtggaC
    """

    def sticky_end_algorithm(x, y, _l):
        return sticky_end_sub_strings(x, y, allow_partial_overlap)

    if allow_blunt:
        algorithm_fn = combine_algorithms(sticky_end_algorithm, blunt_overlap)
    else:
        algorithm_fn = sticky_end_algorithm

    products = common_function_assembly_products(
        frags, None, algorithm_fn, circular_only
    )
    return _recast_sources(products, LigationSource)


def assembly_is_multi_site(asm: list[EdgeRepresentationAssembly]) -> bool:
    """Returns True if the assembly is a multi-site assembly, False otherwise."""

    if len(asm) < 2:
        return False

    is_cycle = asm[0][1] == asm[-1][0]
    asm2 = edge_representation2subfragment_representation(asm, is_cycle)

    return all(f[1] != f[2] for f in asm2)


def gateway_assembly(
    frags: list[Dseqrecord],
    reaction_type: Literal["BP", "LR"],
    greedy: bool = False,
    circular_only: bool = False,
    multi_site_only: bool = False,
) -> list[Dseqrecord]:
    """Returns the products for Gateway assembly / Gateway cloning.

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

    >>> from pydna.assembly2 import gateway_assembly
    >>> from pydna.dseqrecord import Dseqrecord
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


    """

    if reaction_type not in ["BP", "LR"]:
        raise ValueError(
            f"Invalid reaction type: {reaction_type}, can only be BP or LR"
        )

    def algorithm_fn(x, y, _l):
        return gateway_overlap(x, y, reaction_type, greedy)

    filter_results_function = None if not multi_site_only else assembly_is_multi_site

    products = common_function_assembly_products(
        frags, None, algorithm_fn, circular_only, filter_results_function
    )
    products = _recast_sources(
        products,
        GatewaySource,
        reaction_type=reaction_type,
        greedy=greedy,
    )

    if len(products) == 0:
        # Build a list of all the sites in the fragments
        sites_in_fragments = list()
        for frag in frags:
            sites_in_fragments.append(list(find_gateway_sites(frag, greedy).keys()))
        formatted_strings = [
            f'fragment {i + 1}: {", ".join(sites)}'
            for i, sites in enumerate(sites_in_fragments)
        ]
        raise ValueError(
            f"Inputs are not compatible for {reaction_type} reaction.\n\n"
            + "\n".join(formatted_strings),
        )
    return products


def common_function_integration_products(
    frags: list[Dseqrecord], limit: int | None, algorithm: Callable
) -> list[Dseqrecord]:
    """Common function to avoid code duplication for integration products.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to integrate
    limit : int or None
        Minimum overlap length required, or None if not applicable
    algorithm : Callable
        Function that determines valid overlaps between fragments

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules
    """
    if len(frags) == 1:
        asm = SingleFragmentAssembly(frags, limit, algorithm)
    else:
        asm = Assembly(
            frags, limit, algorithm, use_fragment_order=False, use_all_fragments=True
        )

    if frags[0].circular:
        raise ValueError(
            "Genome must be linear for integration assembly, use in vivo assembly instead"
        )

    # We only want insertions in the genome (first fragment)
    output_assemblies = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]
    return [assemble(frags, a, True) for a in output_assemblies]


def common_handle_insertion_fragments(
    genome: Dseqrecord, inserts: list[Dseqrecord]
) -> list[Dseqrecord]:
    """Common function to handle / validate insertion fragments.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord] or Dseqrecord
        DNA fragment(s) to insert

    Returns
    -------
    list[Dseqrecord]
        List containing genome and insert fragments
    """
    if not isinstance(genome, Dseqrecord):
        raise ValueError("Genome must be a Dseqrecord object")

    if not isinstance(inserts, list) or not all(
        isinstance(f, Dseqrecord) for f in inserts
    ):
        raise ValueError("Inserts must be a list of Dseqrecord objects")

    if len(inserts) == 0:
        raise ValueError("Inserts must be a non-empty list of Dseqrecord objects")

    return [genome] + inserts


def common_function_excision_products(
    genome: Dseqrecord, limit: int | None, algorithm: Callable
) -> list[Dseqrecord]:
    """Common function to avoid code duplication for excision products.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    limit : int or None
        Minimum overlap length required, or None if not applicable
    algorithm : Callable
        Function that determines valid overlaps between fragments

    Returns
    -------
    list[Dseqrecord]
        List of excised DNA molecules
    """
    asm = SingleFragmentAssembly([genome], limit, algorithm)
    return asm.assemble_circular() + asm.assemble_insertion()


def homologous_recombination_integration(
    genome: Dseqrecord,
    inserts: list[Dseqrecord],
    limit: int = 40,
) -> list[Dseqrecord]:
    """Returns the products resulting from the integration of an insert (or inserts joined
    through in vivo recombination) into the genome through homologous recombination.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord]
        DNA fragment(s) to insert
    limit : int, optional
        Minimum homology length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules


    Examples
    --------

    Below an example with a single insert.

    >>> from pydna.assembly2 import homologous_recombination_integration
    >>> from pydna.dseqrecord import Dseqrecord
    >>> homology = "AAGTCCGTTCGTTTTACCTG"
    >>> genome = Dseqrecord(f"aaaaaa{homology}ccccc{homology}aaaaaa")
    >>> insert = Dseqrecord(f"{homology}gggg{homology}")
    >>> products = homologous_recombination_integration(genome, [insert], 20)
    >>> str(products[0].seq)
    'aaaaaaAAGTCCGTTCGTTTTACCTGggggAAGTCCGTTCGTTTTACCTGaaaaaa'

    Below an example with two inserts joined through homology.

    >>> homology2 = "ATTACAGCATGGGAAGAAAGA"
    >>> insert_1 = Dseqrecord(f"{homology}gggg{homology2}")
    >>> insert_2 = Dseqrecord(f"{homology2}cccc{homology}")
    >>> products = homologous_recombination_integration(genome, [insert_1, insert_2], 20)
    >>> str(products[0].seq)
    'aaaaaaAAGTCCGTTCGTTTTACCTGggggATTACAGCATGGGAAGAAAGAccccAAGTCCGTTCGTTTTACCTGaaaaaa'
    """
    fragments = common_handle_insertion_fragments(genome, inserts)

    products = common_function_integration_products(
        fragments, limit, common_sub_strings
    )
    return _recast_sources(products, HomologousRecombinationSource)


def homologous_recombination_excision(
    genome: Dseqrecord, limit: int = 40
) -> list[Dseqrecord]:
    """Returns the products resulting from the excision of a fragment from the genome through
    homologous recombination.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    limit : int, optional
        Minimum homology length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List containing excised plasmid and remaining genome sequence

    Examples
    --------

    Example of a homologous recombination event, where a plasmid is excised from the
    genome (circular sequence of 25 bp), and that part is removed from the genome,
    leaving a shorter linear sequence (32 bp).

    >>> from pydna.assembly2 import homologous_recombination_excision
    >>> from pydna.dseqrecord import Dseqrecord
    >>> homology = "AAGTCCGTTCGTTTTACCTG"
    >>> genome = Dseqrecord(f"aaaaaa{homology}ccccc{homology}aaaaaa")
    >>> products = homologous_recombination_excision(genome, 20)
    >>> products
    [Dseqrecord(o25), Dseqrecord(-32)]
    """
    products = common_function_excision_products(genome, limit, common_sub_strings)
    return _recast_sources(products, HomologousRecombinationSource)


def cre_lox_integration(
    genome: Dseqrecord, inserts: list[Dseqrecord]
) -> list[Dseqrecord]:
    """Returns the products resulting from the integration of an insert (or inserts joined
    through cre-lox recombination among them) into the genome through cre-lox integration.

    Also works with lox66 and lox71 (see ``pydna.cre_lox`` for more details).

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord] or Dseqrecord
        DNA fragment(s) to insert

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules

    Examples
    --------

    Below an example of reversible integration and excision.

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import cre_lox_integration, cre_lox_excision
    >>> from pydna.cre_lox import LOXP_SEQUENCE
    >>> a = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
    >>> b = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)
    >>> [a, b]
    [Dseqrecord(-45), Dseqrecord(o39)]
    >>> res = cre_lox_integration(a, [b])
    >>> res
    [Dseqrecord(-84)]
    >>> res2 = cre_lox_excision(res[0])
    >>> res2
    [Dseqrecord(o39), Dseqrecord(-45)]

    Below an example with lox66 and lox71 (irreversible integration).
    Here, the result of excision is still returned because there is a low
    probability of it happening, but it's considered a rare event.

    >>> lox66 = 'ATAACTTCGTATAGCATACATTATACGAACGGTA'
    >>> lox71 = 'TACCGTTCGTATAGCATACATTATACGAAGTTAT'
    >>> a = Dseqrecord(f"cccccc{lox66}aaaaa")
    >>> b = Dseqrecord(f"{lox71}bbbbb", circular=True)
    >>> res = cre_lox_integration(a, [b])
    >>> res
    [Dseqrecord(-84)]
    >>> res2 = cre_lox_excision(res[0])
    >>> res2
    [Dseqrecord(o39), Dseqrecord(-45)]

    """
    fragments = common_handle_insertion_fragments(genome, inserts)
    products = common_function_integration_products(fragments, None, cre_loxP_overlap)
    return _recast_sources(products, CreLoxRecombinationSource)


def cre_lox_excision(genome: Dseqrecord) -> list[Dseqrecord]:
    """Returns the products for CRE-lox excision.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence

    Returns
    -------
    list[Dseqrecord]
        List containing excised plasmid and remaining genome sequence

    Examples
    --------

    Below an example of reversible integration and excision.

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import cre_lox_integration, cre_lox_excision
    >>> from pydna.cre_lox import LOXP_SEQUENCE
    >>> a = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
    >>> b = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)
    >>> [a, b]
    [Dseqrecord(-45), Dseqrecord(o39)]
    >>> res = cre_lox_integration(a, [b])
    >>> res
    [Dseqrecord(-84)]
    >>> res2 = cre_lox_excision(res[0])
    >>> res2
    [Dseqrecord(o39), Dseqrecord(-45)]

    Below an example with lox66 and lox71 (irreversible integration).
    Here, the result of excision is still returned because there is a low
    probability of it happening, but it's considered a rare event.

    >>> lox66 = 'ATAACTTCGTATAGCATACATTATACGAACGGTA'
    >>> lox71 = 'TACCGTTCGTATAGCATACATTATACGAAGTTAT'
    >>> a = Dseqrecord(f"cccccc{lox66}aaaaa")
    >>> b = Dseqrecord(f"{lox71}bbbbb", circular=True)
    >>> res = cre_lox_integration(a, [b])
    >>> res
    [Dseqrecord(-84)]
    >>> res2 = cre_lox_excision(res[0])
    >>> res2
    [Dseqrecord(o39), Dseqrecord(-45)]
    """
    products = common_function_excision_products(genome, None, cre_loxP_overlap)
    return _recast_sources(products, CreLoxRecombinationSource)


def crispr_integration(
    genome: Dseqrecord,
    inserts: list[Dseqrecord],
    guides: list[Primer],
    limit: int = 40,
) -> list[Dseqrecord]:
    """
    Returns the products for CRISPR integration.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord]
        DNA fragment(s) to insert
    guides : list[Primer]
        List of guide RNAs as Primer objects. This may change in the future.
    limit : int, optional
        Minimum overlap length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules

    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly2 import crispr_integration
    >>> from pydna.primer import Primer
    >>> genome = Dseqrecord("aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac", name="genome")
    >>> insert = Dseqrecord("aaccggttAAAAAAAAAttcaaagcac", name="insert")
    >>> guide = Primer("ttcaatgcaaacagtaatga", name="guide")
    >>> product, *_ = crispr_integration(genome, [insert], [guide], 8)
    >>> product
    Dseqrecord(-27)

    """
    if len(guides) == 0:
        raise ValueError("At least one guide RNA is required for CRISPR integration")

    # Get all the possible products from the homologous recombination integration
    products = homologous_recombination_integration(genome, inserts, limit)

    # Verify that the guides cut in the region that will be repaired

    # First we collect the positions where the guides cut
    guide_cuts = []
    for guide in guides:
        enzyme = cas9(str(guide.seq))
        possible_cuts = genome.seq.get_cutsites(enzyme)
        if len(possible_cuts) == 0:
            raise ValueError(
                f"Could not find Cas9 cutsite in the target sequence using the guide: {guide.name}"
            )
        # Keep only the position of the cut
        possible_cuts = [cut[0] for (cut, _) in possible_cuts]
        guide_cuts.append(possible_cuts)

    # Then, we check it the possible homologous recombination products contain the cuts
    # from the guides inside the repair region.
    # We also add the used guides to each product. This is very important!
    valid_products = []
    for i, product in enumerate(products):
        # The second element of product.source.input is conventionally the insert/repair fragment
        # The other two (first and third) are the two bits of the genome
        repair_start = location_boundaries(product.source.input[0].right_location)[0]
        # Here we do +1 because the position of the cut marks the boundary (e.g. 0:10, 10:20 if a cut is at pos 10)
        repair_end = location_boundaries(product.source.input[2].left_location)[1] + 1
        repair_location = create_location(repair_start, repair_end, len(genome))
        some_cuts_inside_repair = []
        all_cuts_inside_repair = []
        for cut_group in guide_cuts:
            cuts_in_repair = [cut for cut in cut_group if cut in repair_location]
            some_cuts_inside_repair.append(len(cuts_in_repair) != 0)
            all_cuts_inside_repair.append(len(cuts_in_repair) == len(cut_group))

        if all(some_cuts_inside_repair):
            used_guides = [g for i, g in enumerate(guides) if all_cuts_inside_repair[i]]
            # Add the used guides to the product <----- VERY IMPORTANT!
            product.source.input.extend([SourceInput(sequence=g) for g in used_guides])
            valid_products.append(product)

            if not all(all_cuts_inside_repair):
                raise ValueError(
                    "Some guides cut outside the repair region, please check the guides"
                )

    if len(valid_products) != len(products):
        warnings.warn(
            "Some recombination products were discarded because they had off-target cuts",
            category=UserWarning,
            stacklevel=2,
        )

    return _recast_sources(valid_products, CRISPRSource)


def pcr_assembly(
    template: Dseqrecord,
    fwd_primer: Primer,
    rvs_primer: Primer,
    add_primer_features: bool = False,
    limit: int = 14,
    mismatches: int = 0,
) -> list[Dseqrecord]:
    """Returns the products for PCR assembly.

    Parameters
    ----------
    template : Dseqrecord
        Template sequence
    fwd_primer : Primer
        Forward primer
    rvs_primer : Primer
        Reverse primer
    add_primer_features : bool, optional
        If True, add primer features to the product, by default False
    limit : int, optional
        Minimum overlap length required, by default 14
    mismatches : int, optional
        Maximum number of mismatches, by default 0

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules
    """

    minimal_annealing = limit + mismatches
    fragments = [fwd_primer, template, rvs_primer]
    asm = PCRAssembly(
        fragments,
        limit=minimal_annealing,
        mismatches=mismatches,
    )
    products = asm.assemble_linear()
    # If both primers are the same, remove duplicates
    if str(fwd_primer.seq).upper() == str(rvs_primer.seq).upper():
        products = [p for p in products if not p.source.input[1].reverse_complemented]
    if add_primer_features:
        products = [annotate_primer_binding_sites(prod, fragments) for prod in products]

    return _recast_sources(products, PCRSource, add_primer_features=add_primer_features)
