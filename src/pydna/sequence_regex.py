# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
import re
from Bio.Data.IUPACData import ambiguous_dna_values

custom_ambiguous_only_dna_values = {**ambiguous_dna_values}
for normal_base in "ACGT":
    del custom_ambiguous_only_dna_values[normal_base]


def compute_regex_site(site: str) -> str:
    """
    Creates a regex pattern from a string that may contain degenerate bases.

    Args:
        site: The string to convert to a regex pattern.

    Returns:
        The regex pattern.
    """
    upper_site = site.upper()
    for k, v in custom_ambiguous_only_dna_values.items():
        if len(v) > 1:
            upper_site = upper_site.replace(k, f"[{''.join(v)}]")

    # Make case insensitive
    upper_site = f"(?i){upper_site}"
    return upper_site


def dseqrecord_finditer(pattern: str, seq: Dseqrecord) -> list[re.Match]:
    """
    Finds all matches of a regex pattern in a Dseqrecord.

    Args:
        pattern: The regex pattern to search for.
        seq: The Dseqrecord to search in.

    Returns:
        A list of matches.
    """
    query = str(seq.seq) if not seq.circular else str(seq.seq) * 2
    matches = re.finditer(pattern, query)
    return (m for m in matches if m.start() <= len(seq))
