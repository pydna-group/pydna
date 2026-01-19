#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module contain functions for primer design for various purposes.

- :func:primer_design for designing primers for a sequence or a matching primer for an existing primer. Returns an :class:`Amplicon` object (same as the :mod:`amplify` module returns).

- :func:assembly_fragments Adds tails to primers for a linear assembly through homologous recombination or Gibson assembly.

- :func:circular_assembly_fragments Adds tails to primers for a circular assembly through homologous recombination or Gibson assembly.

"""

from pydna.tm import tm_default
import math
import copy
from pydna.amplicon import Amplicon
from pydna.amplify import Anneal
from pydna.amplify import pcr
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
import operator
from typing import Tuple
from itertools import pairwise, product
import re


def _design_primer(
    target_tm: float,
    template: Dseqrecord,
    limit: int,
    tm_func,
    starting_length: int = 0,
) -> Tuple[float, str]:
    """returns a tuple (temp, primer)"""

    if starting_length < limit:
        starting_length = limit

    length = starting_length
    tlen = len(template)

    tmp = max(0, tm_func(str(template.seq[:length])))

    p = str(template.seq[:length])

    if tmp < target_tm:
        condition = operator.le
        increment = 1
    else:
        condition = operator.ge
        increment = -1
    while condition(tmp, target_tm):
        prev_temp = tmp
        prev_primer = p
        length += increment
        p = str(template.seq[:length])
        tmp = tm_func(p)
        if length >= tlen:
            break
        # Never go below the limit
        if length < limit:
            return template.seq[:limit]

    if abs(target_tm - tmp) < abs(target_tm - prev_temp):
        return p
    else:
        return prev_primer


def primer_design(
    template,
    fp=None,
    rp=None,
    limit=13,
    target_tm=55.0,
    tm_func=tm_default,
    estimate_function=None,
    **kwargs,
):
    """This function designs a forward primer and a reverse primer for PCR amplification
    of a given template sequence.

    The template argument is a Dseqrecord object or equivalent containing the template sequence.

    The optional fp and rp arguments can contain an existing primer for the sequence (either the forward or reverse primer).
    One or the other primers can be specified, not both (since then there is nothing to design!, use the pydna.amplify.pcr function instead).

    The limit argument is the minimum length of the primer. The default value is 13.

    If one of the primers is given, the other primer is designed to match in terms of Tm.
    If both primers are designed, they will be designed to target_tm

    tm_func is a function that takes an ascii string representing an oligonuceotide as argument and returns a float.
    Some useful functions can be found in the :mod:`pydna.tm` module, but can be substituted for a custom made function.

    estimate_function is a tm_func-like function that is used to get a first guess for the primer design, that is then used as starting
    point for the final result. This is useful when the tm_func function is slow to calculate (e.g. it relies on an
    external API, such as the NEB primer design API). The estimate_function should be faster than the tm_func function.
    The default value is `None`.
    To use the default `tm_func` as estimate function to get the NEB Tm faster, you can do:
    `primer_design(dseqr, target_tm=55, tm_func=tm_neb, estimate_function=tm_default)`.

    The function returns a pydna.amplicon.Amplicon class instance. This object has
    the object.forward_primer and object.reverse_primer properties which contain the designed primers.


    Parameters
    ----------

    template : pydna.dseqrecord.Dseqrecord
        a Dseqrecord object. The only required argument.

    fp, rp : pydna.primer.Primer, optional
        optional pydna.primer.Primer objects containing one primer each.

    target_tm : float, optional
        target tm for the primers, set to 55°C by default.

    tm_func : function
        Function used for tm calculation. This function takes an ascii string
        representing an oligonuceotide as argument and returns a float.
        Some useful functions can be found in the :mod:`pydna.tm` module, but can be
        substituted for a custom made function.

    Returns
    -------
    result : Amplicon

    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> t=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> t
    Dseqrecord(-64)
    >>> from pydna.design import primer_design
    >>> ampl = primer_design(t)
    >>> ampl
    Amplicon(64)
    >>> ampl.forward_primer
    f64 17-mer:5'-atgactgctaacccttc-3'
    >>> ampl.reverse_primer
    r64 18-mer:5'-catcgtaagtttcgaacg-3'
    >>> print(ampl.figure())
    5atgactgctaacccttc...cgttcgaaacttacgatg3
                         ||||||||||||||||||
                        3gcaagctttgaatgctac5
    5atgactgctaacccttc3
     |||||||||||||||||
    3tactgacgattgggaag...gcaagctttgaatgctac5
    >>> pf = "GGATCC" + ampl.forward_primer
    >>> pr = "GGATCC" + ampl.reverse_primer
    >>> pf
    f64 23-mer:5'-GGATCCatgactgct..ttc-3'
    >>> pr
    r64 24-mer:5'-GGATCCcatcgtaag..acg-3'
    >>> from pydna.amplify import pcr
    >>> pcr_prod = pcr(pf, pr, t)
    >>> print(pcr_prod.figure())
          5atgactgctaacccttc...cgttcgaaacttacgatg3
                               ||||||||||||||||||
                              3gcaagctttgaatgctacCCTAGG5
    5GGATCCatgactgctaacccttc3
           |||||||||||||||||
          3tactgacgattgggaag...gcaagctttgaatgctac5
    >>> print(pcr_prod.seq)
    GGATCCatgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatgGGATCC
    >>> from pydna.primer import Primer
    >>> pf = Primer("atgactgctaacccttccttggtgttg", id="myprimer")
    >>> ampl = primer_design(t, fp = pf)
    >>> ampl.forward_primer
    myprimer 27-mer:5'-atgactgctaaccct..ttg-3'
    >>> ampl.reverse_primer
    r64 32-mer:5'-catcgtaagtttcga..atc-3'
    """

    def design(target_tm, template):
        if estimate_function:
            first_guess = _design_primer(target_tm, template, limit, estimate_function)
            return _design_primer(target_tm, template, limit, tm_func, len(first_guess))
        else:
            return _design_primer(target_tm, template, limit, tm_func)

    if not fp and not rp:

        fp = Primer((design(target_tm, template)))
        target_tm = tm_func(str(fp.seq))

        rp = Primer(design(target_tm, template.reverse_complement()))
    elif fp and not rp:
        try:
            fp = Anneal((fp,), template).forward_primers.pop()
        except IndexError:
            raise ValueError("Forward primer does not anneal")
        except Exception:  # pragma: no cover
            print("Unexpected error")
            raise
        target_tm = tm_func(fp.footprint)

        rp = Primer(design(target_tm, template.reverse_complement()))
    elif not fp and rp:
        try:
            rp = Anneal((rp,), template).reverse_primers.pop()
        except IndexError:
            raise ValueError("Reverse primer does not anneal")
        except Exception:  # pragma: no cover
            print("Unexpected error")
            raise
        target_tm = tm_func(rp.footprint)

        fp = Primer(design(target_tm, template))
    else:
        raise ValueError("Specify maximum one of the two primers.")

    if fp.id == "id":  # <unknown id>
        fp.id = "f{}".format(len(template))

    if rp.id == "id":
        rp.id = "r{}".format(len(template))

    if fp.name == "name":
        fp.name = "f{}".format(len(template))

    if rp.name == "name":
        rp.name = "r{}".format(len(template))

    fp.description = fp.id + " " + template.accession
    rp.description = rp.id + " " + template.accession

    ampl = Anneal((fp, rp), template, limit=limit)

    prod = ampl.products[0] if ampl.products else Amplicon("")

    if len(ampl.products) > 1:
        import warnings
        from pydna import _PydnaWarning

        warnings.warn(
            "designed primers do not yield a unique PCR product", _PydnaWarning
        )

    return prod


def assembly_fragments(f, overlap=35, maxlink=40, circular=False):
    """This function return a list of :mod:`pydna.amplicon.Amplicon` objects where
    primers have been modified with tails so that the fragments can be fused in
    the order they appear in the list by for example Gibson assembly or homologous
    recombination.

    Given that we have two linear :mod:`pydna.amplicon.Amplicon` objects a and b

    we can modify the reverse primer of a and forward primer of b with tails to allow
    fusion by fusion PCR, Gibson assembly or in-vivo homologous recombination.
    The basic requirements for the primers for the three techniques are the same.

    ::



            _________ a _________
           /                     \\
           agcctatcatcttggtctctgca
                             |||||
                            <gacgt
           agcct>
           |||||
           tcggatagtagaaccagagacgt

                                   __________ b ________
                                  /                     \\
                                  TTTATATCGCATGACTCTTCTTT
                                                    |||||
                                                   <AGAAA
                                  TTTAT>
                                  |||||
                                  AAATATAGCGTACTGAGAAGAAA


           agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
           ||||||||||||||||||||||||||||||||||||||||||||||
           tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
           \\___________________ c ______________________/


    Design tailed primers incorporating a part of the next or previous fragment to be assembled.

    ::


      agcctatcatcttggtctctgca
      |||||||||||||||||||||||
                      gagacgtAAATATA

      |||||||||||||||||||||||
      tcggatagtagaaccagagacgt


                             TTTATATCGCATGACTCTTCTTT
                             |||||||||||||||||||||||

                      ctctgcaTTTATAT
                             |||||||||||||||||||||||
                             AAATATAGCGTACTGAGAAGAAA

    PCR products with flanking sequences are formed in the PCR process.

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATA
                      \\____________/

                         identical
                         sequences
                       ____________
                      /            \\
                      ctctgcaTTTATATCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    The fragments can be fused by any of the techniques mentioned earlier to form c:

    ::

      agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
      ||||||||||||||||||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA


    The first argument of this function is a list of sequence objects containing
    Amplicons and other similar objects.

    **At least every second sequence object needs to be an Amplicon**

    This rule exists because if a sequence object is that is not a PCR product
    is to be fused with another fragment, that other fragment needs to be an Amplicon
    so that the primer of the other object can be modified to include the whole stretch
    of sequence homology needed for the fusion. See the example below where a is a
    non-amplicon (a linear plasmid  vector for instance)

    ::

       _________ a _________           __________ b ________
      /                     \\         /                     \\
      agcctatcatcttggtctctgca   <-->  TTTATATCGCATGACTCTTCTTT
      |||||||||||||||||||||||         |||||||||||||||||||||||
      tcggatagtagaaccagagacgt                          <AGAAA
                                      TTTAT>
                                      |||||||||||||||||||||||
                                <-->  AAATATAGCGTACTGAGAAGAAA


           agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
           ||||||||||||||||||||||||||||||||||||||||||||||
           tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
           \\___________________ c ______________________/


    In this case only the forward primer of b is fitted with a tail with a part a:

    ::


      agcctatcatcttggtctctgca
      |||||||||||||||||||||||
      tcggatagtagaaccagagacgt


                             TTTATATCGCATGACTCTTCTTT
                             |||||||||||||||||||||||
                                              <AGAAA
               tcttggtctctgcaTTTATAT
                             |||||||||||||||||||||||
                             AAATATAGCGTACTGAGAAGAAA

    PCR products with flanking sequences are formed in the PCR process.

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATA
                      \\____________/

                         identical
                         sequences
                       ____________
                      /            \\
                      ctctgcaTTTATATCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    The fragments can be fused by for example Gibson assembly:

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaacca

                                   TCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    to form c:

    ::

      agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
      ||||||||||||||||||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA


    The first argument of this function is a list of sequence objects containing
    Amplicons and other similar objects.

    The overlap argument controls how many base pairs of overlap required between
    adjacent sequence fragments. In the junction between Amplicons, tails with the
    length of about half of this value is added to the two primers
    closest to the junction.

    ::

            >       <
            Amplicon1
                     Amplicon2
                     >       <

                     ⇣

            >       <-
            Amplicon1
                     Amplicon2
                    ->       <

    In the case of an Amplicon adjacent to a Dseqrecord object, the tail will
    be twice as long (1*overlap) since the
    recombining sequence is present entirely on this primer:

    ::

            Dseqrecd1
                     Amplicon1
                     >       <

                     ⇣

            Dseqrecd1
                     Amplicon1
                   -->       <

    Note that if the sequence of DNA fragments starts or stops with an Amplicon,
    the very first and very last prinmer will not be modified i.e. assembles are
    always assumed to be linear. There are simple tricks around that for circular
    assemblies depicted in the last two examples below.

    The maxlink arguments controls the cut off length for sequences that will be
    synhtesized by adding them to primers for the adjacent fragment(s). The
    argument list may contain short spacers (such as spacers between fusion proteins).


    ::

        Example 1: Linear assembly of PCR products (pydna.amplicon.Amplicon class objects) ------


        >       <         >       <
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon4
                 >       <         >       <

                             ⇣
                             pydna.design.assembly_fragments
                             ⇣

        >       <-       ->       <-                      pydna.assembly.Assembly
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon4     ➤  Amplicon1Amplicon2Amplicon3Amplicon4
                ->       <-       ->       <


        Example 2: Linear assembly of alternating Amplicons and other fragments


        >       <         >       <
        Amplicon1         Amplicon2
                 Dseqrecd1         Dseqrecd2

                             ⇣
                             pydna.design.assembly_fragments
                             ⇣

        >       <--     -->       <--                     pydna.assembly.Assembly
        Amplicon1         Amplicon2
                 Dseqrecd1         Dseqrecd2     ➤  Amplicon1Dseqrecd1Amplicon2Dseqrecd2


        Example 3: Linear assembly of alternating Amplicons and other fragments


        Dseqrecd1         Dseqrecd2
                 Amplicon1         Amplicon2
                 >       <       -->       <

                             ⇣
                     pydna.design.assembly_fragments
                             ⇣
                                                          pydna.assembly.Assembly
        Dseqrecd1         Dseqrecd2
                 Amplicon1         Amplicon2     ➤  Dseqrecd1Amplicon1Dseqrecd2Amplicon2
               -->       <--     -->       <


        Example 4: Circular assembly of alternating Amplicons and other fragments

                         ->       <==
        Dseqrecd1         Amplicon2
                 Amplicon1         Dseqrecd1
               -->       <-
                             ⇣
                             pydna.design.assembly_fragments
                             ⇣
                                                           pydna.assembly.Assembly
                         ->       <==
        Dseqrecd1         Amplicon2                    -Dseqrecd1Amplicon1Amplicon2-
                 Amplicon1                       ➤    |                             |
               -->       <-                            -----------------------------

        ------ Example 5: Circular assembly of Amplicons

        >       <         >       <
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon1
                 >       <         >       <

                             ⇣
                             pydna.design.assembly_fragments
                             ⇣

        >       <=       ->       <-
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon1
                ->       <-       +>       <

                             ⇣
                     make new Amplicon using the Amplicon1.template and
                     the last fwd primer and the first rev primer.
                             ⇣
                                                           pydna.assembly.Assembly
        +>       <=       ->       <-
         Amplicon1         Amplicon3                  -Amplicon1Amplicon2Amplicon3-
                  Amplicon2                      ➤   |                             |
                 ->       <-                          -----------------------------




    Parameters
    ----------

    f : list of :mod:`pydna.amplicon.Amplicon` and other Dseqrecord like objects
        list Amplicon and Dseqrecord object for which fusion primers should be constructed.

    overlap : int, optional
        Length of required overlap between fragments.

    maxlink : int, optional
        Maximum length of spacer sequences that may be present in f. These will be included in tails for designed primers.

    circular : bool, optional
        If True, the assembly is circular. If False, the assembly is linear.

    Returns
    -------
    seqs : list of :mod:`pydna.amplicon.Amplicon` and other Dseqrecord like objects :mod:`pydna.amplicon.Amplicon` objects

        ::

          [Amplicon1,
           Amplicon2, ...]


    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.design import primer_design
    >>> a=primer_design(Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg"))
    >>> b=primer_design(Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc"))
    >>> c=primer_design(Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt"))
    >>> from pydna.design import assembly_fragments
    >>> # We would like a circular recombination, so the first sequence has to be repeated
    >>> fa1,fb,fc,fa2 = assembly_fragments([a,b,c,a])
    >>> # Since all fragments are Amplicons, we need to extract the rp of the 1st and fp of the last fragments.
    >>> from pydna.amplify import pcr
    >>> fa = pcr(fa2.forward_primer, fa1.reverse_primer, a)
    >>> [fa,fb,fc]
    [Amplicon(100), Amplicon(101), Amplicon(102)]
    >>> fa.name, fb.name, fc.name = "fa fb fc".split()
    >>> from pydna.assembly import Assembly
    >>> assemblyobj = Assembly([fa,fb,fc])
    >>> assemblyobj
    Assembly
    fragments....: 100bp 101bp 102bp
    limit(bp)....: 25
    G.nodes......: 6
    algorithm....: common_sub_strings
    >>> assemblyobj.assemble_linear()
    [Contig(-231), Contig(-166), Contig(-36)]
    >>> assemblyobj.assemble_circular()[0].seguid()
    'cdseguid=85t6tfcvWav0wnXEIb-lkUtrl4s'
    >>> (a+b+c).looped().seguid()
    'cdseguid=85t6tfcvWav0wnXEIb-lkUtrl4s'
    >>> print(assemblyobj.assemble_circular()[0].figure())
     -|fa|36
    |     \\/
    |     /\\
    |     36|fb|36
    |           \\/
    |           /\\
    |           36|fc|36
    |                 \\/
    |                 /\\
    |                 36-
    |                    |
     --------------------
    >>>

    """

    # Recursive call for circular assemblies
    if circular:
        fragments = assembly_fragments(
            f + f[0:1], overlap=overlap, maxlink=maxlink, circular=False
        )

        if hasattr(fragments[0], "template"):
            fragments[0] = pcr(
                (fragments[-1].forward_primer, fragments[0].reverse_primer),
                fragments[0].template,
            )
        return fragments[:-1]

    # sanity check for arguments
    nf = [item for item in f if len(item) > maxlink]
    if not all(
        hasattr(i[0], "template") or hasattr(i[1], "template") for i in zip(nf, nf[1:])
    ):
        raise ValueError(
            "Every second fragment larger than maxlink has to be an Amplicon object"
        )

    f = [copy.copy(f) for f in f]

    first_fragment_length = len(f[0])
    last_fragment_length = len(f[-1])

    if first_fragment_length <= maxlink:
        # first fragment should be removed and added to second fragment (new first fragment) forward primer
        f[1].forward_primer = f[0].seq._data.decode("ASCII") + f[1].forward_primer

        f = f[1:]
    # else:

    if last_fragment_length <= maxlink:
        f[-2].reverse_primer = (
            f[-1].seq.reverse_complement()._data.decode("ASCII") + f[-2].reverse_primer
        )
        f = f[:-1]

    # else:

    empty = Dseqrecord("")

    tail_length = math.ceil(overlap / 2)

    for i in range(len(f) - 1):
        first = f[i]
        secnd = f[i + 1]

        secnd_len = len(secnd)

        if secnd_len <= maxlink:

            third = f[i + 2]

            if hasattr(f[i], "template") and hasattr(third, "template"):

                #    "secnd is is flanked by amplicons, so half of secnd should be added each flanking primers"
                # )

                first.reverse_primer = (
                    secnd.seq.reverse_complement()._data.decode("ASCII")[
                        secnd_len // 2 :
                    ]
                    + first.reverse_primer
                )
                third.forward_primer = (
                    secnd.seq._data.decode("ASCII")[secnd_len // 2 :]
                    + third.forward_primer
                )

                lnk = (
                    third.seq.reverse_complement()._data.decode("ASCII")
                    + secnd.reverse_complement().seq._data.decode("ASCII")[
                        : secnd_len // 2
                    ]
                )[-tail_length:]

                first.reverse_primer = lnk + first.reverse_primer

                lnk = (
                    first.seq._data.decode("ASCII")
                    + secnd.seq._data.decode("ASCII")[: secnd_len // 2]
                )[-tail_length:]

                third.forward_primer = lnk + third.forward_primer

            elif hasattr(first, "template"):
                first.reverse_primer = (
                    secnd.seq.reverse_complement()._data.decode("ASCII")
                    + first.reverse_primer
                )
                lnk = str(third.seq[:overlap].reverse_complement())
                first.reverse_primer = lnk + first.reverse_primer
            elif hasattr(third, "template"):
                third.forward_primer = (
                    secnd.seq._data.decode("ASCII") + third.forward_primer
                )
                lnk = str(first.seq[-overlap:])
                third.forward_primer = lnk + third.forward_primer
            secnd = empty
            f[i + 2] = third
        else:  # secnd is larger than maxlink
            if hasattr(first, "template") and hasattr(secnd, "template"):
                lnk = str(first.seq[-tail_length:])
                # #_module_logger.debug("4 %s", lnk)
                secnd.forward_primer = lnk + secnd.forward_primer
                lnk = str(secnd.seq[:tail_length].reverse_complement())
                # #_module_logger.debug("5 %s", lnk)
                first.reverse_primer = lnk + first.reverse_primer
            elif hasattr(first, "template"):
                lnk = str(secnd.seq[:overlap].reverse_complement())
                # #_module_logger.debug("6 %s", lnk)
                first.reverse_primer = lnk + first.reverse_primer
            elif hasattr(secnd, "template"):
                lnk = str(first.seq[-overlap:])
                # #_module_logger.debug("7 %s", lnk)
                secnd.forward_primer = lnk + secnd.forward_primer
        f[i] = first
        f[i + 1] = secnd

    f = [item for item in f if len(item)]

    return [
        (
            pcr(
                p.forward_primer,
                p.reverse_primer,
                p.template,
                limit=min((p.forward_primer._fp, p.reverse_primer._fp)),
            )
            if hasattr(p, "template")
            else p
        )
        for p in f
    ]


def circular_assembly_fragments(f, overlap=35, maxlink=40):
    """
    Equivalent to `assembly_fragments` with `circular=True`.

    Deprecated, kept for backward compatibility. Use `assembly_fragments` with `circular=True` instead.
    """
    import warnings

    warnings.warn(
        "The circular_assembly_fragments function is deprecated. Use assembly_fragments with circular=True instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return assembly_fragments(f, overlap=overlap, maxlink=maxlink, circular=True)


def user_assembly_design(
    f: list[Amplicon], max_overlap: int = 15, min_overlap: int = 4, max_tail=50
) -> list[Amplicon]:

    import warnings

    warnings.warn(
        "The user_assembly_design function is experimental and "
        "may change in future versions.",
        category=FutureWarning,
        stacklevel=2,
    )

    assert max_overlap > min_overlap, (
        f"max_overlap ({max_overlap}) "
        "has to be larger than min_overlap "
        f"({min_overlap})"
    )
    amplicons = []

    for fragment in f:
        amplicons.append(primer_design(fragment))

    flag = True

    for ths, nxt in pairwise(amplicons):

        A_positions_in_ths = [m.start() for m in re.finditer("A|a", str(ths.seq))]
        T_positions_in_nxt = [m.start() for m in re.finditer("T|t", str(nxt.seq))]

        for ths_a, ths_t in zip(A_positions_in_ths[::-1], T_positions_in_nxt):

            sticky_length = ths_t + len(ths) - ths_a

            if sticky_length < min_overlap:
                continue

            if sticky_length > max_overlap:
                flag = False
                break

            rp = bytearray(
                nxt.seq[: ths_t + 1].rc()._data + ths.reverse_primer.seq._data
            )
            fp = bytearray(ths.seq[ths_a:]._data + nxt.forward_primer.seq._data)

            fp[sticky_length] = ord(b"U")
            rp[sticky_length] = ord(b"U")

            ths.reverse_primer = Primer(rp)
            nxt.forward_primer = Primer(fp)

            break  # Primers were designed.
        else:
            flag = False

        if flag:
            continue

        # No suitable T-A pair was found on opposite sides of both fragments
        # Look for T-A pairs contained in either sequence
        # Distance between the T-A and proximity to the junction are important
        # factors
        T_positions_in_ths = [m.start() for m in re.finditer("T|t", str(ths.seq))]
        pairs_ths = product(A_positions_in_ths[::-1], T_positions_in_ths[::-1])

        for ths_a, ths_t in pairs_ths:
            if ths_a > ths_t:
                continue
            sticky_length = ths_t - ths_a
            if sticky_length < min_overlap:
                continue
            if sticky_length > max_overlap:
                continue
            pair_ths = ths_a, ths_t
            break
        else:
            pair_ths = tuple()
            ths_a, ths_t = 0, 0

        A_positions_in_nxt = [m.start() for m in re.finditer("A|a", str(nxt.seq))]
        pairs_nxt = product(A_positions_in_nxt, T_positions_in_nxt)

        for nxt_a, nxt_t in pairs_nxt:
            if nxt_a > nxt_t:
                continue
            sticky_length = nxt_t - nxt_a
            if sticky_length < min_overlap:
                continue
            if sticky_length > max_overlap:
                continue
            pair_nxt = nxt_a, nxt_t
            break
        else:
            pair_nxt = tuple()
            nxt_a, nxt_t = 0, 0

        if (pair_ths and not pair_nxt) or len(ths) - ths_a <= nxt_t:
            # T-A pair in ths;
            # Move ths reverse primer downstream
            # Extend nxt foward primer tail

            fp = bytearray(ths.seq[ths_a:]._data + nxt.forward_primer.seq._data)
            fp[ths_t - ths_a] = ord(b"U")
            nxt.forward_primer = Primer(fp)
            shorter_ths = ths[: ths_t + 1]
            rp = bytearray(
                primer_design(
                    shorter_ths, limit=ths_t - ths_a + 1
                ).reverse_primer.seq._data
            )
            rp[ths_t - ths_a] = ord(b"U")
            ths.reverse_primer = Primer(rp)

        elif (not pair_ths and pair_nxt) or len(ths) - ths_a >= nxt_t:
            # T-A pair in nxt; modify ths reverse primer
            # Move nxt forward primer upstream
            # Extend ths reverse primer tail
            rp = bytearray(
                nxt.seq[: nxt_t + 1].rc()._data + ths.reverse_primer.seq._data
            )
            rp[nxt_t - nxt_a] = ord(b"U")
            ths.reverse_primer = Primer(rp)
            shorter_nxt = nxt[nxt_a:]
            fp = bytearray(
                primer_design(
                    shorter_nxt, limit=nxt_t - nxt_a + 1
                ).forward_primer.seq._data
            )
            fp[nxt_t - nxt_a] = ord(b"U")
            nxt.forward_primer = Primer(fp)

    return amplicons
