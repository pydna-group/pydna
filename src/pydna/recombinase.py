# -*- coding: utf-8 -*-
"""
Generic site-specific recombinase functionality.

This module provides tools for simulating site-specific recombination reactions
mediated by recombinases (e.g. phage integrases such as phiC31, Bxb1, etc.).

Recombinase sites are specified as strings where the **lowercase** portion
indicates the homology/overlap core shared between the two sites. The
uppercase portion represents the flanking recognition arms.

For example::

    >>> site1 = "ATGCCCTAAaaTT"
    >>> site2 = "AAaaTTTTTTTCCCT"

The lowercase ``aa`` is the overlap that will appear in the assembled product.

Sites may contain IUPAC degenerate bases (N, W, S, etc.) and the search
handles both linear and circular sequences.

The easiest way to use it is with the recombinase_integration and
recombinase_excision functions from the assembly2 module.

Integration and excision with a single recombinase::

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.recombinase import Recombinase
    >>> from pydna.assembly2 import recombinase_integration, recombinase_excision
    >>> site1 = "ATGCCCTAAaaCT"
    >>> site2 = "CAaaTTTTTTTCCCT"
    >>> genome = Dseqrecord("ccccccATGCCCTAAAACTaaaaa")
    >>> insert = Dseqrecord("CAAATTTTTTTCCCTbbbbb", circular=True)
    >>> rec = Recombinase(site1, site2)
    >>> products = recombinase_integration(genome, [insert], rec)
    >>>
    >>> # Cre-Lox style (same site on both sides): integrate then excise returns originals
    >>> site = "ATGaaGTA"
    >>> genome = Dseqrecord("ccccccATGAAGTAAAAA")
    >>> insert = Dseqrecord("ATGAAGTABbbbb", circular=True)
    >>> rec = Recombinase(site, site)
    >>> integrated = recombinase_integration(genome, [insert], rec)[0]
    >>> excised = recombinase_excision(integrated, rec)

Find and annotate sites in a sequence::

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.recombinase import Recombinase
    >>> site1, site2 = "ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT"
    >>> rec = Recombinase(site1, site2, site1_name="mysite1", site2_name="mysite2")
    >>> seq = Dseqrecord("gggATGCCCTAAaaTTttt")
    >>> sites = rec.find(seq)
    >>> annotated = rec.annotate(seq)

When several sites are possible, or when using multiple recombinases, you can
create a RecombinaseCollection. It has the methods overlap, find, and annotate,
so you can use it as a single Recombinase. For an example, check the gateway module::

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.recombinase import Recombinase, RecombinaseCollection
    >>> rec1 = Recombinase("AAaaTTC", "CCaaTTC", site1_name="s1", site2_name="s2")
    >>> rec2 = Recombinase("GAccACC", "TCccAAC", site1_name="s3", site2_name="s4")
    >>> collection = RecombinaseCollection([rec1, rec2])
    >>> seq = Dseqrecord("gggAAAATTCTtttGACCACCTttt")
    >>> sites = collection.find(seq)
    >>> annotated = collection.annotate(seq)

IUPAC degenerate bases (e.g. N matches any base)::

    >>> site1 = "ATGNNNaaTT"
    >>> site2 = "CCNNaaTTGG"
    >>> rec = Recombinase(site1, site2)

Using a Recombinase as Assembly algorithm::

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.recombinase import Recombinase
    >>> from pydna.assembly2 import Assembly
    >>> site1 = "ATGCCCTAAaaTT"
    >>> site2 = "AAaaTTTTTTTCCCT"
    >>> seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    >>> seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")
    >>> rec = Recombinase(site1, site2)
    >>> asm = Assembly([seqA, seqB], algorithm=rec.overlap, use_fragment_order=False)
    >>> products = asm.assemble_linear()

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
    >>> _ = rec.overlap(seqA, seqB)
    >>> sites = rec.find(seqA)
    >>> _ = rec.annotate(seqA)
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
            (self._site2_fwd, self._site1_fwd, off2, off1),
            (self._site2_rev, self._site1_rev, off2_rev, off1_rev),
        ]

    def __repr__(self) -> str:
        return f"Recombinase(site1={self.site1}, site2={self.site2}, site1_name={self.site1_name}, site2_name={self.site2_name})"

    def get_reverse_recombinase(
        self, site12_name: str = "site12", site21_name: str = "site21"
    ) -> "Recombinase":
        """Get a recombinase that does the reverse reaction.

        Parameters
        ----------
        site12_name : str, optional
            Label for site12 in find/annotate output. Default is "site12".
        site21_name : str, optional
            Label for site21 in find/annotate output. Default is "site21".

        Returns
        -------
        Recombinase
            A recombinase that does the reverse reaction.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> from pydna.recombinase import Recombinase
        >>> rec = Recombinase("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
        >>> rec.get_reverse_recombinase()
        Recombinase(site1=ATGCCCTAAaaTTTTTTTCCCT, site2=AAaaTT, site1_name=site12, site2_name=site21)
        """
        _, _, off1, off2 = self._configs[0]
        site12 = (
            self.site1[: off1 + self._homology_len]
            + self.site2[off2 + self._homology_len :]
        )
        site21 = (
            self.site2[: off2 + self._homology_len]
            + self.site1[off1 + self._homology_len :]
        )
        return Recombinase(site12, site21, site12_name, site21_name)

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


class RecombinaseCollection:
    """A collection of recombinases."""

    def __init__(self, recombinases: list[Recombinase]):
        if not isinstance(recombinases, list):
            raise ValueError("recombinases must be a list of Recombinase objects")
        if not all(isinstance(r, Recombinase) for r in recombinases):
            raise ValueError("recombinases must be a list of Recombinase objects")
        if len(recombinases) == 0:
            raise ValueError("recombinases must be a non-empty list")
        self.recombinases = recombinases

    def overlap(self, seqx: Dseqrecord, seqy: Dseqrecord) -> list[SequenceOverlap]:
        """Find overlaps between *seqx* and *seqy* mediated by the recombinases."""
        return sum((r.overlap(seqx, seqy) for r in self.recombinases), [])

    def find(self, seq: Dseqrecord) -> dict[str, list[SimpleLocation]]:
        """Find all occurrences of the recombinase sites in *seq*."""
        out = dict()
        for rec in self.recombinases:
            rec_sites = rec.find(seq)
            for k, v in rec_sites.items():
                out.setdefault(k, []).extend(v)
        return out

    def annotate(self, seq: Dseqrecord) -> Dseqrecord:
        """Annotate *seq* with features for all occurrences of the recombinase sites."""
        out = seq
        for rec in self.recombinases:
            out = rec.annotate(out)
        return out
