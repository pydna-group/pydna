# -*- coding: utf-8 -*-
"""
This module contains the functions for oligonucleotide hybridization.
"""

from pydna.common_sub_strings import common_sub_strings
from Bio.Seq import reverse_complement
from pydna.primer import Primer
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.opencloning_models import OligoHybridizationSource, SourceInput


def oligonucleotide_hybridization_overhangs(
    fwd_oligo_seq: str, rvs_oligo_seq: str, minimal_annealing: int
) -> list[int]:
    """
    Returns possible overhangs between two oligos given a minimal annealing length, and
    returns an error if mismatches are found.

    see https://github.com/manulera/OpenCloning_backend/issues/302 for notation

    >>> from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization_overhangs
    >>> oligonucleotide_hybridization_overhangs("ATGGC", "GCCAT", 3)
    [0]
    >>> oligonucleotide_hybridization_overhangs("aATGGC", "GCCAT", 5)
    [-1]
    >>> oligonucleotide_hybridization_overhangs("ATGGC", "GCCATa", 5)
    [1]
    >>> oligonucleotide_hybridization_overhangs("ATGGC", "GCCATaaGCCAT", 5)
    [0, 7]

    If the minimal annealing length is longer than the length of the shortest oligo, it returns an empty list.

    >>> oligonucleotide_hybridization_overhangs("ATGGC", "GCCATaaGCCAT", 100)
    []

    If it's possible to anneal for ``minimal_annealing`` length, but with mismatches, it raises an error.

    >>> oligonucleotide_hybridization_overhangs("cATGGC", "GCCATa", 5)
    Traceback (most recent call last):
        ...
    ValueError: The oligonucleotides can anneal with mismatches
    """
    matches = common_sub_strings(
        fwd_oligo_seq.lower(),
        reverse_complement(rvs_oligo_seq.lower()),
        minimal_annealing,
    )

    for pos_fwd, pos_rvs, length in matches:

        if (pos_fwd != 0 and pos_rvs != 0) or (
            pos_fwd + length < len(fwd_oligo_seq)
            and pos_rvs + length < len(rvs_oligo_seq)
        ):
            raise ValueError("The oligonucleotides can anneal with mismatches")

    # Return possible overhangs
    return [pos_rvs - pos_fwd for pos_fwd, pos_rvs, length in matches]


def oligonucleotide_hybridization(
    fwd_primer: Primer, rvs_primer: Primer, minimal_annealing: int
) -> list[Dseqrecord]:
    """
    Returns a list of Dseqrecord objects representing the hybridization of two primers.

    >>> from pydna.primer import Primer
    >>> from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization
    >>> fwd_primer = Primer("ATGGC")
    >>> rvs_primer = Primer("GCCA")
    >>> oligonucleotide_hybridization(fwd_primer, rvs_primer, 3)[0].seq
    Dseq(-5)
    ATGGC
     ACCG

    Multiple values can be returned:

    >>> rvs_primer2 = Primer("GCCATaaGCCAT")
    >>> oligonucleotide_hybridization(fwd_primer, rvs_primer2, 3)[0].seq
    Dseq(-12)
    ATGGC
    TACCGaaTACCG
    >>> oligonucleotide_hybridization(fwd_primer, rvs_primer2, 3)[1].seq
    Dseq(-12)
           ATGGC
    TACCGaaTACCG

    If no possible overhangs are found, it returns an empty list.

    >>> oligonucleotide_hybridization(fwd_primer, rvs_primer, 100)
    []

    If there are mismatches given the minimal annealing length, it raises an error.

    >>> fwd_primer3 = Primer("cATGGC")
    >>> rvs_primer3 = Primer("GCCATa")
    >>> oligonucleotide_hybridization(fwd_primer3, rvs_primer3, 5)
    Traceback (most recent call last):
        ...
    ValueError: The oligonucleotides can anneal with mismatches
    """
    possible_overhangs = oligonucleotide_hybridization_overhangs(
        str(fwd_primer.seq), str(rvs_primer.seq), minimal_annealing
    )
    sources = [
        OligoHybridizationSource(
            overhang_crick_3prime=pos,
            input=[SourceInput(sequence=fwd_primer), SourceInput(sequence=rvs_primer)],
        )
        for pos in possible_overhangs
    ]
    return [
        Dseqrecord(
            Dseq(
                str(fwd_primer.seq),
                str(rvs_primer.seq),
                ovhg=source.overhang_crick_3prime,
            ),
            source=source,
        )
        for source in sources
    ]
