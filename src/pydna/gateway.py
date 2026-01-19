# -*- coding: utf-8 -*-
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
import re
import itertools
from Bio.SeqFeature import SimpleLocation, SeqFeature
from pydna.utils import shift_location
from pydna.sequence_regex import compute_regex_site, dseqrecord_finditer


raw_gateway_common = {
    "attB1": "CHWVTWTGTACAAAAAANNNG",
    "attB2": "CHWVTWTGTACAAGAAANNNG",
    "attB3": "CHWVTWTGTATAATAAANNNG",
    "attB4": "CHWVTWTGTATAGAAAANNNG",
    "attB5": "CHWVTWTGTATACAAAANNNG",
    "attL1": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAAAAANNNG",
    "attL2": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAGAAANNNG",
    "attL3": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAATAAANNNG",
    "attL4": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAGAAAANNNG",
    "attL5": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATACAAAANNNG",
    "attR1": "CHWVTWTGTACAAAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attR2": "CHWVTWTGTACAAGAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attR3": "CHWVTWTGTATAATAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attR4": "CHWVTWTGTATAGAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attR5": "CHWVTWTGTATACAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "overlap_1": "twtGTACAAAaaa",
    "overlap_2": "twtGTACAAGaaa",
    "overlap_3": "twtGTATAATaaa",
    "overlap_4": "twtGTATAGAaaa",
    "overlap_5": "twtGTATACAaaa",
}


raw_gateway_sites_greedy = {
    **raw_gateway_common,
    "attP1": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attP2": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAGAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attP3": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAATAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attP4": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAGAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
    "attP5": "VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATACAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV",
}

raw_gateway_sites_conservative = {
    **raw_gateway_common,
    "attP1": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAAAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP2": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAGAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP3": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTATAATAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP4": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTATAGAAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP5": "AAAWWAWKRWTTTWWTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTATACAAAAGYWGAACGAGAARCGTAAARTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAAYCAACTACTTAGATGGTATTAGTGACCTGTA",
}

gateway_sites_greedy = {
    k: {
        "forward_regex": compute_regex_site(v),
        "reverse_regex": compute_regex_site(reverse_complement(v)),
        "consensus_sequence": v,
    }
    for k, v in raw_gateway_sites_greedy.items()
}

gateway_sites_conservative = {
    k: {
        "forward_regex": compute_regex_site(v),
        "reverse_regex": compute_regex_site(reverse_complement(v)),
        "consensus_sequence": v,
    }
    for k, v in raw_gateway_sites_conservative.items()
}

# From snapgene - ask Valerie
primer_design_attB = {
    "attB1": "ACAAGTTTGTACAAAAAAGCAGGCT",
    "attB2": "ACCACTTTGTACAAGAAAGCTGGGT",
    "attB3": "ACAACTTTGTATAATAAAGTTGTA",
    "attB4": "ACAACTTTGTATAGAAAAGTTGTA",
    "attB5": "ACAACTTTGTATACAAAAGTTGTA",
}


def gateway_overlap(
    seqx: Dseqrecord, seqy: Dseqrecord, reaction: str, greedy: bool
) -> list[tuple[int, int, int]]:
    """
    Find gateway overlaps. If greedy is True, it uses a more greedy consensus site to find attP sites,
    which might give false positives
    """
    if reaction not in ["BP", "LR"]:
        raise ValueError(f"Invalid overlap type: {reaction}")

    gateway_sites = gateway_sites_greedy if greedy else gateway_sites_conservative
    out = list()
    # Iterate over the four possible att sites
    for num in range(1, 5):
        # Iterate over the two possible orientations
        # The sites have to be in the same orientation (fwd + fwd or rev + rev)
        for pattern in ["forward_regex", "reverse_regex"]:
            # The overlap regex is the same for all types
            overlap_regex = gateway_sites[f"overlap_{num}"][pattern]

            # Iterate over pairs B, P and P, B for BP and L, R and R, L for LR
            for site_x, site_y in zip(reaction, reaction[::-1]):
                site_x_regex = gateway_sites[f"att{site_x}{num}"][pattern]
                matches_x = list(dseqrecord_finditer(site_x_regex, seqx))
                if len(matches_x) == 0:
                    continue

                site_y_regex = gateway_sites[f"att{site_y}{num}"][pattern]
                matches_y = list(dseqrecord_finditer(site_y_regex, seqy))
                if len(matches_y) == 0:
                    continue

                for match_x, match_y in itertools.product(matches_x, matches_y):
                    # Find the overlap sequence within each match, and use the
                    # core 7 pbs that are constant
                    overlap_x = re.search(overlap_regex, match_x.group())
                    overlap_y = re.search(overlap_regex, match_y.group())

                    # Sanity check
                    assert (
                        overlap_x is not None and overlap_y is not None
                    ), "Something went wrong, no overlap found within the matches"

                    out.append(
                        (
                            match_x.start() + overlap_x.start() + 3,
                            match_y.start() + overlap_y.start() + 3,
                            7,
                        )
                    )

    return out


def find_gateway_sites(
    seq: Dseqrecord, greedy: bool
) -> dict[str, list[SimpleLocation]]:
    """Find all gateway sites in a sequence and return a dictionary with the name and positions of the sites."""
    gateway_sites = gateway_sites_greedy if greedy else gateway_sites_conservative
    out = dict()
    for site in gateway_sites:
        if not site.startswith("att"):
            continue

        for pattern in ["forward_regex", "reverse_regex"]:
            matches = list(dseqrecord_finditer(gateway_sites[site][pattern], seq))
            for match in matches:
                if site not in out:
                    out[site] = []
                strand = 1 if pattern == "forward_regex" else -1
                loc = SimpleLocation(match.start(), match.end(), strand)
                loc = shift_location(loc, 0, len(seq))
                out[site].append(loc)
    return out


def annotate_gateway_sites(seq: Dseqrecord, greedy: bool) -> Dseqrecord:
    sites = find_gateway_sites(seq, greedy)
    for site in sites:
        for loc in sites[site]:
            seq.features.append(
                SeqFeature(loc, type="protein_bind", qualifiers={"label": [site]})
            )
    return seq
