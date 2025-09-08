# -*- coding: utf-8 -*-
from itertools import product
from pydna.dseqrecord import Dseqrecord
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import reverse_complement
from pydna.sequence_regex import compute_regex_site, dseqrecord_finditer
from Bio.SeqFeature import Location, SimpleLocation, SeqFeature
from pydna.utils import shift_location

# We create a dictionary to map ambiguous bases to their consensus base
# For example, ambigous_base_dict['ACGT'] -> 'N'
ambiguous_base_dict = {}
for ambiguous, bases in ambiguous_dna_values.items():
    ambiguous_base_dict["".join(sorted(bases))] = ambiguous

# To handle N values
ambiguous_base_dict["N"] = "N"

# This is the original loxP sequence, here for reference
LOXP_SEQUENCE = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"

loxP_sequences = [
    # https://blog.addgene.org/plasmids-101-cre-lox
    # loxP
    "ATAACTTCGTATANNNTANNNTATACGAAGTTAT",
    # PMID:12202778
    # lox66
    "ATAACTTCGTATANNNTANNNTATACGAACGGTA",
    # lox71
    "TACCGTTCGTATANNNTANNNTATACGAAGTTAT",
]

loxP_consensus = ""

for pos in range(len(LOXP_SEQUENCE)):
    all_letters = set(seq[pos] for seq in loxP_sequences)
    key = "".join(sorted(all_letters))
    loxP_consensus += ambiguous_base_dict[key]

# We compute the regex for the forward and reverse loxP sequences
loxP_regex = (
    compute_regex_site(loxP_consensus),
    compute_regex_site(reverse_complement(loxP_consensus)),
)


def cre_loxP_overlap(
    x: Dseqrecord, y: Dseqrecord, _l: None = None
) -> list[tuple[int, int, int]]:
    """Find matching loxP sites between two sequences."""
    out = list()
    for pattern in loxP_regex:
        matches_x = dseqrecord_finditer(pattern, x)
        matches_y = dseqrecord_finditer(pattern, y)

        for match_x, match_y in product(matches_x, matches_y):
            value_x = match_x.group()
            value_y = match_y.group()
            if value_x[13:21] == value_y[13:21]:
                out.append((match_x.start() + 13, match_y.start() + 13, 8))
    # Unique values (keeping the order)
    unique_out = []
    for item in out:
        if item not in unique_out:
            unique_out.append(item)
    return unique_out


loxP_dict = {
    "loxP": "ATAACTTCGTATANNNTANNNTATACGAAGTTAT",
    "lox66": "ATAACTTCGTATANNNTANNNTATACGAACGGTA",
    "lox71": "TACCGTTCGTATANNNTANNNTATACGAAGTTAT",
    "loxP_mutant": "TACCGTTCGTATANNNTANNNTATACGAACGGTA",
}


def get_regex_dict(original_dict: dict[str, str]) -> dict[str, str]:
    """Get the regex dictionary for the original dictionary."""
    out = dict()
    for site in original_dict:
        consensus_seq = original_dict[site]
        is_palindromic = consensus_seq == reverse_complement(consensus_seq)
        out[site] = {
            "forward_regex": compute_regex_site(original_dict[site]),
            "reverse_regex": (
                None
                if is_palindromic
                else compute_regex_site(reverse_complement(original_dict[site]))
            ),
        }
    return out


def find_loxP_sites(seq: Dseqrecord) -> dict[str, list[Location]]:
    """Find all loxP sites in a sequence and return a dictionary with the name and positions of the sites."""

    out = dict()
    regex_dict = get_regex_dict(loxP_dict)
    for site in loxP_dict:

        for pattern in ["forward_regex", "reverse_regex"]:
            # Palindromic sequences have no reverse complement
            if regex_dict[site][pattern] is None:
                continue
            matches = list(dseqrecord_finditer(regex_dict[site][pattern], seq))
            for match in matches:
                if site not in out:
                    out[site] = []
                strand = 1 if pattern == "forward_regex" else -1
                loc = SimpleLocation(match.start(), match.end(), strand)
                loc = shift_location(loc, 0, len(seq))
                out[site].append(loc)
    return out


def annotate_loxP_sites(seq: Dseqrecord) -> Dseqrecord:
    sites = find_loxP_sites(seq)
    for site in sites:
        for loc in sites[site]:
            # Don't add the same feature twice
            if not any(
                f.location == loc
                and f.type == "protein_bind"
                and f.qualifiers.get("label", []) == [site]
                for f in seq.features
            ):
                seq.features.append(
                    SeqFeature(loc, type="protein_bind", qualifiers={"label": [site]})
                )
    return seq
