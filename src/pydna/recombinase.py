# -*- coding: utf-8 -*-
"""
Generic site-specific recombinase functionality.

This module provides tools for simulating site-specific recombination reactions
mediated by recombinases (e.g. phage integrases such as phiC31, Bxb1, etc.).

Recombinase sites are specified as strings where the **lowercase** portion
indicates the homology/overlap core shared between the two sites. The
uppercase portion represents the flanking recognition arms.

For example::

    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

The lowercase ``aa`` is the overlap that will appear in the assembled product.

Sites may contain IUPAC degenerate bases (N, W, S, etc.) and the search
handles both linear and circular sequences.
"""

import re
import itertools
import copy

from Bio.Seq import reverse_complement
from Bio.SeqFeature import SimpleLocation, SeqFeature

from pydna.dseqrecord import Dseqrecord
from pydna.utils import shift_location
from pydna.sequence_regex import compute_regex_site, dseqrecord_finditer
from pydna.types import SequenceOverlap


def _recombinase_homology_offset_and_length(site: str) -> tuple[int, int]:
    """
    Return (offset, length) of the lowercase homology region inside a
    recombinase recognition site string.

    The lowercase segment represents the shared homology between the two
    recombinase sites.

    Parameters
    ----------
    site : str
        The recognition site sequence, with the homology region in lowercase.
        Must contain only letters (no digits or symbols).

    Returns
    -------
    offset : int
        Index of the first lowercase character in the site.
    length : int
        Length of the contiguous lowercase homology region.

    Raises
    ------
    ValueError
        If the site contains invalid characters or no lowercase region.
    """
    if not re.fullmatch(r"[A-Z]+[a-z]+[A-Z]+", site):
        raise ValueError(
            "Recombinase recognition site is not in the expected format."
            "Expected format: [A-Z]+[a-z]+[A-Z]+, e.g. 'ATGCCCTAAaaTT'"
        )
    match = re.search(r"[a-z]+", site)
    return match.start(), len(match.group())


class Recombinase:
    """A site-specific recombinase defined by two recognition sites. Homology cores must match.

    Parameters
    ----------
    site1 : str
        The first recognition site. The homology core must be in lowercase.
    site2 : str
        The second recognition site. The homology core must be in lowercase.
    site1_name : str, optional
        Label for site1 in find/annotate output. Default is "site1".
    site2_name : str, optional
        Label for site2 in find/annotate output. Default is "site2".

    Examples
    --------
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.recombinase import Recombinase
    >>> rec = Recombinase("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
    >>> seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    >>> seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")
    >>> rec.overlap(seqA, seqB)
    [(12, 6, 2)]
    >>> sites = rec.find(seqA)
    >>> rec.annotate(seqA)
    """

    def __init__(
        self,
        site1: str,
        site2: str,
        site1_name: str = "site1",
        site2_name: str = "site2",
    ):
        self.site1 = site1
        self.site2 = site2
        self.site1_name = site1_name
        self.site2_name = site2_name

        off1, len1 = _recombinase_homology_offset_and_length(site1)
        off2, len2 = _recombinase_homology_offset_and_length(site2)
        if site1[off1 : off1 + len1] != site2[off2 : off2 + len2]:
            raise ValueError(
                "Recombinase recognition sites do not have matching homology cores."
                f"Expected {site1[off1:off1 + len1]} == {site2[off2:off2 + len2]}"
            )
        self._homology_len = len1

        rc1 = reverse_complement(site1)
        rc2 = reverse_complement(site2)
        off1_rev, _ = _recombinase_homology_offset_and_length(rc1)
        off2_rev, _ = _recombinase_homology_offset_and_length(rc2)

        self._site1_fwd = compute_regex_site(site1)
        self._site1_rev = compute_regex_site(rc1)
        self._site2_fwd = compute_regex_site(site2)
        self._site2_rev = compute_regex_site(rc2)

        self._configs = [
            (self._site1_fwd, self._site2_fwd, off1, off2),
            (self._site1_rev, self._site2_rev, off1_rev, off2_rev),
        ]

    def overlap(
        self,
        seqx: Dseqrecord,
        seqy: Dseqrecord,
        limit: None = None,
    ) -> list[SequenceOverlap]:
        """Find overlaps between *seqx* and *seqy* mediated by the sites.

        Parameters
        ----------
        seqx : Dseqrecord
            The first sequence.
        seqy : Dseqrecord
            The second sequence.
        limit : None
            Not used. Present only to satisfy the
            ``AssemblyAlgorithmType`` convention.

        Returns
        -------
        list[SequenceOverlap]
            A list of ``(start_in_x, start_in_y, length)`` tuples.
        """
        matches: list[SequenceOverlap] = []

        for site_x_re, site_y_re, off_x, off_y in self._configs:
            matches_x = list(dseqrecord_finditer(site_x_re, seqx))
            if not matches_x:
                continue
            matches_y = list(dseqrecord_finditer(site_y_re, seqy))
            if not matches_y:
                continue

            for mx, my in itertools.product(matches_x, matches_y):
                # Here we wrap in case the sequence is circular and the match spans the origin,
                # and the offset makes the match go beyond the origin.
                matches.append(
                    (
                        (mx.start() + off_x) % len(seqx),
                        (my.start() + off_y) % len(seqy),
                        self._homology_len,
                    )
                )

        # Deduplicate while preserving order
        seen: set[SequenceOverlap] = set()
        unique: list[SequenceOverlap] = []
        for m in matches:
            if m not in seen:
                seen.add(m)
                unique.append(m)
        return unique

    def find(self, seq: Dseqrecord) -> dict[str, list[SimpleLocation]]:
        """Find all occurrences of the recombinase sites in *seq*.

        Parameters
        ----------
        seq : Dseqrecord
            Sequence to search.

        Returns
        -------
        dict[str, list[SimpleLocation]]
            Dictionary mapping site names to lists of locations.
        """
        out: dict[str, list[SimpleLocation]] = {}
        for name, raw_site in [
            (self.site1_name, self.site1),
            (self.site2_name, self.site2),
        ]:
            fwd_re = compute_regex_site(raw_site)
            rev_re = compute_regex_site(reverse_complement(raw_site))
            for pattern, strand in [(fwd_re, 1), (rev_re, -1)]:
                for match in dseqrecord_finditer(pattern, seq):
                    loc = SimpleLocation(match.start(), match.end(), strand)
                    loc = shift_location(loc, 0, len(seq))
                    out.setdefault(name, []).append(loc)
        return out

    def annotate(self, seq: Dseqrecord) -> Dseqrecord:
        """Annotate *seq* with features for all occurrences of the recombinase sites.

        Parameters
        ----------
        seq : Dseqrecord
            Sequence to annotate (modified in place and returned).

        Returns
        -------
        Dseqrecord
            The same sequence with added features.
        """
        out_seq = copy.deepcopy(seq)
        sites = self.find(out_seq)
        for name, locs in sites.items():
            for loc in locs:
                if not any(
                    f.location == loc
                    and f.type == "protein_bind"
                    and f.qualifiers.get("label", []) == [name]
                    for f in out_seq.features
                ):
                    out_seq.features.append(
                        SeqFeature(
                            loc, type="protein_bind", qualifiers={"label": [name]}
                        )
                    )
        return out_seq
