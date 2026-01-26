#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Provides the Dseq class for handling double stranded DNA sequences.

Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class
which can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology.
"""

import itertools
import re
import copy
import sys
import math
import inspect
from typing import List, Tuple, Union

from Bio.Restriction import RestrictionBatch
from Bio.Restriction import CommOnly

from seguid import ldseguid
from seguid import cdseguid

from pydna.seq import Seq
from Bio.Seq import _SeqAbstractBaseClass
from Bio.Data.IUPACData import unambiguous_dna_weights
from Bio.Data.IUPACData import unambiguous_rna_weights
from Bio.Data.IUPACData import atom_weights
from pydna._pretty import pretty_str
from pydna.utils import rc
from pydna.utils import flatten
from pydna.utils import cuts_overlap

from pydna.alphabet import basepair_dict
from pydna.alphabet import dscode_to_watson_table
from pydna.alphabet import dscode_to_crick_table
from pydna.alphabet import regex_ds_melt_factory
from pydna.alphabet import regex_ss_melt_factory
from pydna.alphabet import dscode_to_full_sequence_table
from pydna.alphabet import dscode_to_watson_tail_table
from pydna.alphabet import dscode_to_crick_tail_table
from pydna.alphabet import complement_table_for_dscode
from pydna.alphabet import letters_not_in_dscode
from pydna.alphabet import get_parts
from pydna.alphabet import representation_tuple
from pydna.alphabet import dsbreaks

from pydna.common_sub_strings import common_sub_strings
from pydna.types import DseqType, EnzymesType, CutSiteType


# Sequences larger than this gets a truncated representation.
length_limit_for_repr = 30
placeholder = letters_not_in_dscode[-1]


class CircularBytes(bytes):
    """
    A circular bytes sequence: indexing and slicing wrap around index 0.
    """

    def __new__(cls, value: bytes | bytearray | memoryview):
        return super().__new__(cls, bytes(value))

    def __getitem__(self, key):
        n = len(self)
        if n == 0:
            if isinstance(key, slice):
                return self.__class__(b"")
            raise IndexError("CircularBytes index out of range (empty bytes)")

        if isinstance(key, int):
            return super().__getitem__(key % n)

        if isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            step = 1 if step is None else step
            if step == 0:
                raise ValueError("slice step cannot be zero")

            if step > 0:
                start = 0 if start is None else start
                stop = n if stop is None else stop
                while stop <= start:
                    stop += n
                rng = range(start, stop, step)
            else:
                start = (n - 1) if start is None else start
                stop = -1 if stop is None else stop
                while stop >= start:
                    stop -= n
                rng = range(start, stop, step)

            limit = n if step % n == 0 else n * 2
            out = bytearray()
            count = 0
            for i in rng:
                out.append(super().__getitem__(i % n))
                count += 1
                if count > limit:
                    break
            return self.__class__(bytes(out))

        return super().__getitem__(key)

    def cutaround(self, start: int, length: int) -> bytes:
        """
        Return a circular slice of given length starting at index `start`.
        Can exceed len(self), wrapping around as needed.

        Examples
        --------
        s = CircularBytes(b"ABCDE")
        assert s.cutaround(3, 7) == b"DEABCDE"
        assert s.cutaround(-1, 4) == b"EABC"
        """
        n = len(self)
        if n == 0 or length <= 0:
            return self.__class__(b"")

        start %= n
        out = bytearray()
        for i in range(length):
            out.append(self[(start + i) % n])
        return self.__class__(bytes(out))

    def find(
        self,
        sub: bytes | bytearray | memoryview | str,
        start: int = 0,
        end: int | None = None,
    ) -> int:
        """
        Find a subsequence in the circular sequence, possibly
        wrapping across the origin.
        Returns -1 if not found.
        """
        n = len(self)
        if n == 0:
            return -1

        end = n if end is None else min(end, n)
        doubled = self + self
        try:
            sub = sub.encode("ascii")
        except AttributeError:
            pass

        pos = doubled.find(bytes(sub), start, n + len(sub) - 1)

        if pos == -1 or pos >= n:
            return -1
        return pos


class Dseq(Seq):
    """Dseq describes a double stranded DNA fragment, linear or circular.

    Dseq can be initiated in two ways, using two strings, each representing the
    Watson (upper, sense) strand, the Crick (lower, antisense) strand and an
    optional value describing the stagger betwen the strands on the left side (ovhg).

    Alternatively, a single string represenation using dsIUPAC codes can be used.
    If a single string is used, the letters of that string are interpreted as base
    pairs rather than single bases. For example "A" would indicate the basepair
    "A/T". An expanded IUPAC code is used where the letters PEXI have been assigned
    to GATC on the Watson strand with no paring base on the Crick strand G/"", A/"",
    T/"" and C/"". The letters QFZJ have been assigned the opposite base pairs with
    an empty Watson strand ""/G, ""/A, ""/T, and ""/C.

    ::

        PEXIGATCQFZJ  would indicate the linear double-stranded fragment:

        GATCGATC
            CTAGCTAG



    Parameters
    ----------
    watson : str
        a string representing the Watson (sense) DNA strand or a basepair
        represenation.

    crick : str, optional
        a string representing the Crick (antisense) DNA strand.

    ovhg : int, optional
        A positive or negative number to describe the stagger between the
        Watson and Crick strands.
        see below for a detailed explanation.

    circular : bool, optional
        True indicates that sequence is circular, False that it is linear.


    Examples
    --------
    Dseq is a subclass of the Biopython Bio.Seq.Seq class. The constructor
    can accept two strings representing the Watson (sense) and Crick(antisense)
    DNA strands. These are interpreted as single stranded DNA. There is a check
    for complementarity between the strands.

    If the DNA molecule is staggered on the left side, an integer ovhg
    (overhang) must be given, describing the stagger between the Watson and Crick strand
    in the 5' end of the fragment.

    Additionally, the optional boolean parameter circular can be given to indicate if the
    DNA molecule is circular.

    The most common usage of the Dseq class is probably not to use it directly, but to
    create it as part of a Dseqrecord object (see :class:`pydna.dseqrecord.Dseqrecord`).
    This works in the same way as for the relationship between the :class:`Bio.Seq.Seq` and
    :class:`Bio.SeqRecord.SeqRecord` classes in Biopython.

    There are multiple ways of creating a Dseq object directly listed below, but you can also
    use the function Dseq.from_full_sequence_and_overhangs() to create a Dseq:

    Two arguments (string, string), no overhang provided:

    >>> from pydna.dseq import Dseq
    >>> Dseq("gggaaat","ttt")
    Dseq(-7)
    gggaaat
       ttt

    If Watson and Crick are given, but not ovhg, an attempt will be made to find the best annealing
    between the strands. There are important limitations to this. If there are several ways to
    anneal the strands, this will fail. For long fragments it is quite slow.

    Three arguments (string, string, ovhg=int):

    The ovhg parameter is an integer describing the length of the Crick strand overhang on the
    left side (the 5' end of Watson strand).

    The ovhg parameter controls the stagger at the five prime end::

        dsDNA       overhang

          nnn...    2
        nnnnn...

         nnnn...    1
        nnnnn...

        nnnnn...    0
        nnnnn...

        nnnnn...   -1
         nnnn...

        nnnnn...   -2
          nnn...

    Example of creating Dseq objects with different amounts of stagger:

    >>> Dseq(watson="att", crick="acata", ovhg=-2)
    Dseq(-7)
    att
      ataca
    >>> Dseq(watson="ata",crick="acata",ovhg=-1)
    Dseq(-6)
    ata
     ataca
    >>> Dseq(watson="taa",crick="actta",ovhg=0)
    Dseq(-5)
    taa
    attca
    >>> Dseq(watson="aag",crick="actta",ovhg=1)
    Dseq(-5)
     aag
    attca
    >>> Dseq(watson="agt",crick="actta",ovhg=2)
    Dseq(-5)
      agt
    attca

    If the ovhg parameter is specified a Crick strand also needs to be supplied, or
    an exception is raised.

    >>> Dseq(watson="agt", ovhg=2)
    Traceback (most recent call last):
        ...
    ValueError: ovhg (overhang) defined without a crick strand.


    The shape or topology of the fragment is set by the circular parameter, True or False (default).

    >>> Dseq("aaa", "ttt", ovhg = 0)  # A linear sequence by default
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa", "ttt", ovhg = 0, circular = False)  # A linear sequence if circular is False
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa", "ttt", ovhg = 0, circular = True)  # A circular sequence
    Dseq(o3)
    aaa
    ttt
    >>> Dseq("aaa", "ttt", ovhg=1, circular = False)
    Dseq(-4)
     aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=-1)
    Dseq(-4)
    aaa
     ttt
    >>> Dseq("aaa", "ttt", circular = True , ovhg=0)
    Dseq(o3)
    aaa
    ttt

    >>> a=Dseq("tttcccc","aaacccc")
    >>> a
    Dseq(-11)
        tttcccc
    ccccaaa
    >>> a.ovhg
    4

    >>> b=Dseq("ccccttt","ccccaaa")
    >>> b
    Dseq(-11)
    ccccttt
        aaacccc
    >>> b.ovhg
    -4
    >>>


    dsIUPAC [#]_ is an nn extension to the IUPAC alphabet used to describe ss regions:

    ::

            aaaGATC       GATCccc          ad-hoc representations
        CTAGttt               gggCTAG

        QFZJaaaPEXI       PEXIcccQFZJ      dsIUPAC



    Coercing to string

    >>> str(a)
    'ggggtttcccc'

    A Dseq object can be longer that either the watson or crick strands.

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG

        <-- length -->
              GATCCTTT
        AAAGCCCTA

    The slicing of a linear Dseq object works mostly as it does for a string.

    >>> s="ggatcc"
    >>> s[2:3]
    'a'
    >>> s[2:4]
    'at'
    >>> s[2:4:-1]
    ''
    >>> s[::2]
    'gac'
    >>> from pydna.dseq import Dseq
    >>> d=Dseq(s, circular=False)
    >>> d[2:3]
    Dseq(-1)
    a
    t
    >>> d[2:4]
    Dseq(-2)
    at
    ta
    >>> d[2:4:-1]
    Dseq(-0)
    <BLANKLINE>
    <BLANKLINE>
    >>> d[::2]
    Dseq(-3)
    gac
    ctg


    The slicing of a circular Dseq object has a slightly different meaning.


    >>> s="ggAtCc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggAtCc
    ccTaGg
    >>> d[4:3]
    Dseq(-5)
    CcggA
    GgccT


    The slice [X:X] produces an empty slice for a string, while this
    will return the linearized sequence starting at X:

    >>> s="ggatcc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggatcc
    cctagg
    >>> d[3:3]
    Dseq(-6)
    tccgga
    aggcct
    >>>


    See Also
    --------
    pydna.dseqrecord.Dseqrecord

    """

    def __init__(
        self,
        watson: Union[str, bytes],
        crick: Union[str, bytes, None] = None,
        ovhg=None,
        circular=False,
        pos=0,
    ):
        if isinstance(watson, (bytes, bytearray)):
            # watson is decoded to a string if needed.
            watson = watson.decode("ascii")
        if isinstance(crick, (bytes, bytearray)):
            # crick is decoded to a string if needed.
            crick = crick.decode("ascii")

        if crick is None:
            if ovhg is not None:
                raise ValueError("ovhg (overhang) defined without a crick strand.")
            """
            Giving only the watson string implies inferring the Crick complementary strand
            from the Watson sequence. The watson string can contain dscode letters wich will
            be interpreted as outlined in the pydna.alphabet module.

            The _data property must be a byte string for compatibility with
            Biopython Bio.Seq.Seq
            """
            data = watson
            self._data = data.encode("ascii")

        else:
            """
            Crick strand given, ovhg is optional. An important consequence is that the
            watson and crick strands are interpreted as single stranded DNA that is
            supposed to anneal.

            If ovhg was not given, we try to guess the value below. This will fail
            if there are two or more ways to anneal with equal length of the double
            stranded part.
            """
            if ovhg is None:  # ovhg not given, try to guess from sequences
                limit = int(math.log(len(watson)) / math.log(4))
                olaps = common_sub_strings(
                    str(watson).lower(),
                    str(rc(crick).lower()),
                    limit,
                )

                """No overlaps found, strands do not anneal"""
                if len(olaps) == 0:
                    raise ValueError(
                        "Could not anneal the two strands."
                        f" looked for annealing with at least {limit} basepairs"
                        " Please provide and overhang value (ovhg parameter)"
                    )

                """
                We extract the positions and length of the first (longest) overlap,
                since common_sub_strings sorts the overlaps by length, longest first.
                """

                (pos_watson, pos_crick, longest_olap_length), *rest = olaps

                """
                We see if there is another overlap of the same length
                This means that annealing is ambigous. User should provide
                and ovhg value.
                """
                if any(
                    olap_length >= longest_olap_length for _, _, olap_length in rest
                ):
                    raise ValueError(
                        "More than one way of annealing the"
                        " strands. Please provide ovhg value"
                    )

                ovhg = pos_crick - pos_watson

            """
            Pad both strands on left side ovhg spaces
            a negative number gives no padding,
            """
            sense = ovhg * " " + watson
            antisense = -ovhg * " " + crick[::-1]

            max_len = max(len(sense), len(antisense))

            """pad both strands on right side to same size."""
            sense = sense.ljust(max_len)
            antisense = antisense.ljust(max_len)
            """both strands padded so that bsepairs align"""
            assert len(sense) == len(antisense)

            data = []

            for w, c in zip(sense, antisense):
                try:
                    data.append(basepair_dict[w, c])
                except KeyError as err:
                    print(f"Base mismatch in representation {err}")
                    raise ValueError(f"Base mismatch in representation: {err}")
            data = "".join(data).strip()
            self._data = data.encode("ascii")

        self.circular = circular
        self.pos = pos

        if circular:
            data += data[0:1]

        dsb = dsbreaks(data)

        if dsb:
            msg = "".join(dsb)
            raise ValueError(
                f"Molecule is internally split in {len(dsb)} location(s):\n\n{msg}".strip()
            )

    @classmethod
    def quick(cls, data: bytes, *args, circular=False, pos=0, **kwargs):
        """Fastest way to instantiate an object of the Dseq class.

        No checks of parameters are made.
        Does not call Bio.Seq.Seq.__init__() which has lots of time consuming checks.
        """
        obj = cls.__new__(cls)
        obj.circular = circular
        obj.pos = pos
        obj._data = data

        return obj

    @classmethod
    def from_representation(cls, dsdna: str, *args, **kwargs):
        obj = cls.__new__(cls)
        obj.circular = False
        obj.pos = 0
        clean = inspect.cleandoc("\n" + dsdna)
        watson, crick = [
            ln
            for ln in clean.splitlines()
            if ln.strip() and not ln.strip().startswith("Dseq(")
        ]
        ovhgw = len(watson) - len(watson.lstrip())
        ovhgc = -(len(crick) - len(crick.lstrip()))

        ovhg = ovhgw or ovhgc

        watson = watson.strip()
        crick = crick.strip()[::-1]

        return Dseq(watson, crick, ovhg)

    @classmethod
    def from_full_sequence_and_overhangs(
        cls, full_sequence: str, crick_ovhg: int, watson_ovhg: int
    ):
        """Create a linear Dseq object from a full sequence and the 3' overhangs of each strand.

        The order of the parameters is like this because the 3' overhang of the crick strand is the one
        on the left side of the sequence.


        Parameters
        ----------
        full_sequence: str
            The full sequence of the Dseq object.

        crick_ovhg: int
            The overhang of the crick strand in the 3' end. Equivalent to Dseq.ovhg.

        watson_ovhg: int
            The overhang of the watson strand in the 5' end.

        Returns
        -------
        Dseq
            A Dseq object.

        Examples
        --------

        >>> Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=2)
        Dseq(-6)
          AAAA
        TTTT
        >>> Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=2)
        Dseq(-6)
        AAAAAA
          TT
        >>> Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=-2)
        Dseq(-6)
          AA
        TTTTTT
        >>> Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=-2)
        Dseq(-6)
        AAAA
          TTTT

        """
        full_sequence_rev = str(Dseq(full_sequence).reverse_complement())
        watson = full_sequence
        crick = full_sequence_rev

        # If necessary, we trim the left side
        if crick_ovhg < 0:
            crick = crick[:crick_ovhg]
        elif crick_ovhg > 0:
            watson = watson[crick_ovhg:]

        # If necessary, we trim the right side
        if watson_ovhg < 0:
            watson = watson[:watson_ovhg]
        elif watson_ovhg > 0:
            crick = crick[watson_ovhg:]

        return Dseq(watson, crick=crick, ovhg=crick_ovhg)

    @property
    def watson(self) -> str:
        """
        The watson (upper) strand of the double stranded fragment 5'-3'.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self._data.decode("ascii").translate(dscode_to_watson_table).strip()

    @property
    def crick(self) -> str:
        """
        The crick (lower) strand of the double stranded fragment 5'-3'.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self._data.decode("ascii").translate(dscode_to_crick_table).strip()[::-1]

    @property
    def left_ovhg(self) -> int:
        """
        The 5' overhang of the lower strand compared the the upper.

        See module docstring for more information.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        parts = self.get_parts()
        if parts.single_watson or parts.single_crick:
            return None
        return -len(parts.sticky_left5) or len(parts.sticky_left3)

    ovhg = left_ovhg

    @property
    def right_ovhg(self) -> int:
        """Overhang at the right side (end)."""
        parts = self.get_parts()
        if parts.single_watson or parts.single_crick:
            return None
        return -len(parts.sticky_right5) or len(parts.sticky_right3)

    watson_ovhg = right_ovhg

    def __str__(self) -> str:
        """
        A string representation of the sequence. The returned string
        is the watson strand of a blunt version of the sequence.

        >>> ds = Dseq.from_representation(
        ... '''
        ... GAATTC
        ...   TAA
        ... ''')

        >>> str(ds)
        'GAATTC'
        >>> ds = Dseq.from_representation(
        ... '''
        ...   ATT
        ... CTTAAG
        ... ''')

        >>> str(ds)
        'GAATTC'

        Returns
        -------
        str
            A string representation of the sequence.

        """
        return bytes(self).decode("ascii")

    to_blunt_string = __str__  # alias of __str__ # TODO: consider removing

    def __bytes__(self) -> bytes:
        return self._data.translate(dscode_to_full_sequence_table)

    def mw(self) -> float:
        """The molecular weight of the DNA/RNA molecule in g/mol.

        The molecular weight data in Biopython Bio.Data.IUPACData
        is used. The DNA is assumed to have a 5'-phosphate as many
        DNA fragments from restriction digestion do:

        ::

             P - G-A-T-T-A-C-A - OH
                 | | | | | | |
            OH - C-T-A-A-T-G-T - P

        The molecular weights listed in the unambiguous_dna_weights
        dictionary refers to free monophosphate nucleotides.
        One water molecule is removed for every phopshodiester bond
        formed between nucleotides. For linear molecules, the weight
        of one water molecule is added to account for the terminal
        hydroxyl group and a hydrogen on the 5' terminal phosphate
        group.

        ::

             P - G---A---T - OH  P - C---A - OH
                 |   |   |           |   |
            OH - C---T---A---A---T---G---T - P

        If the DNA is discontinuous, the internal 5'- end is assumed
        to have a phosphate and the 3'- a hydroxyl group:


        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds_lin_obj = Dseq("GATTACA")
        >>> ds_lin_obj
        Dseq(-7)
        GATTACA
        CTAATGT
        >>> round(ds_lin_obj.mw(), 1)
        4359.8
        >>> ds_circ_obj = Dseq("GATTACA", circular = True)
        >>> round(ds_circ_obj.mw(), 1)
        4323.8
        >>> ssobj = Dseq("PEXXEIE")
        >>> ssobj
        Dseq(-7)
        GATTACA
        <BLANKLINE>
        >>> round(ssobj.mw(), 1)
        2184.4
        >>> ds_lin_obj2 = Dseq("GATZFCA")
        >>> ds_lin_obj2
        Dseq(-7)
        GAT  CA
        CTAATGT
        >>> round(ds_lin_obj2.mw(), 1)
        3724.4
        """

        h2o = atom_weights["H"] * 2 + atom_weights["O"]

        mwd = unambiguous_rna_weights | unambiguous_dna_weights | {" ": 0}

        watsn_weight = sum(mwd[nt] - h2o for nt in self.watson.upper())
        crick_weight = sum(mwd[nt] - h2o for nt in self.crick.upper())

        watsn_weight += h2o * len(re.findall(r" +", self.watson))
        crick_weight += h2o * len(re.findall(r" +", self.crick))

        if watsn_weight and not self.circular:
            watsn_weight += h2o

        if crick_weight and not self.circular:
            crick_weight += h2o

        return watsn_weight + crick_weight

    def find(
        self, sub: Union[_SeqAbstractBaseClass, str, bytes], start=0, end=sys.maxsize
    ) -> int:
        """This method behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Returns -1 if the subsequence is NOT found.

        The search is case sensitive.

        Parameters
        ----------

        sub : string or Seq object
            a string or another Seq object to look for.

        start : int, optional
            slice start.

        end : int, optional
            slice end.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> seq = Dseq("agtaagt")
        >>> seq
        Dseq(-7)
        agtaagt
        tcattca
        >>> seq.find("taa")
        2
        >>> seq = Dseq(watson="agta",crick="actta",ovhg=-2)
        >>> seq
        Dseq(-7)
        agta
          attca
        >>> seq.find("taa")
        -1
        >>> seq = Dseq(watson="agta",crick="actta",ovhg=-2)
        >>> seq
        Dseq(-7)
        agta
          attca
        >>> seq.find("ta")
        2
        """
        if self.circular:
            result = CircularBytes(self._data).find(sub, start, end)
        else:
            result = super().find(sub, start, end)
        return result

    def __contains__(self, sub: [str, bytes]) -> bool:
        return self.find(sub) != -1

    def __getitem__(self, sl: [slice, int]) -> DseqType:
        if isinstance(sl, int):
            sl = slice(sl, sl + 1, 1)
        sl = slice(sl.start, sl.stop, sl.step)
        if self.circular:
            cb = CircularBytes(self._data)
            return self.quick(cb[sl])
        return super().__getitem__(sl)

    def __eq__(self, other: DseqType) -> bool:
        """Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties.

        This comparison is case insensitive.

        """
        try:
            same = (
                other.watson.lower() == self.watson.lower()
                and other.crick.lower() == self.crick.lower()
                and other.ovhg == self.ovhg
                and self.circular == other.circular
            )
            # Also test for alphabet ?
        except AttributeError:
            same = False
        return same

    def __repr__(self, lim: int = length_limit_for_repr) -> pretty_str:

        header = f"{self.__class__.__name__}({({False: '-', True: 'o'}[self.circular])}{len(self)})"

        w, c = representation_tuple(
            self._data.decode("ascii"), length_limit_for_repr=length_limit_for_repr
        )

        return pretty_str(header + "\n" + w + "\n" + c)

    def reverse_complement(self) -> "Dseq":
        """Dseq object where watson and crick have switched places.

        This represents the same double stranded sequence.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc")
        >>> a
        Dseq(-8)
        catcgatc
        gtagctag
        >>> b=a.reverse_complement()
        >>> b
        Dseq(-8)
        gatcgatg
        ctagctac
        >>>

        """
        return Dseq.quick(rc(self._data), circular=self.circular)

    rc = reverse_complement  # alias for reverse_complement

    def shifted(self: DseqType, shift: int) -> DseqType:
        """
        Shifted copy of a circular Dseq object.

        >>> ds = Dseq("TAAG", circular = True)
        >>> ds.shifted(1) # First bp moved to right side:
        Dseq(o4)
        AAGT
        TTCA
        >>> ds.shifted(-1) # Last bp moved to left side:
        Dseq(o4)
        GTAA
        CATT
        """
        if not self.circular:
            raise TypeError("DNA is not circular.")
        shift = shift % len(self)
        if not shift:
            return copy.deepcopy(self)
        else:
            return (self[shift:] + self[:shift]).looped()

    def looped(self: DseqType) -> DseqType:
        """Circularized Dseq object.

        This can only be done if the two ends are compatible,
        otherwise a TypeError is raised.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc")
        >>> a
        Dseq(-8)
        catcgatc
        gtagctag
        >>> a.looped()
        Dseq(o8)
        catcgatc
        gtagctag
        >>> b = Dseq("iatcgatj")
        >>> b
        Dseq(-8)
        catcgat
         tagctag
        >>> b.looped()
        Dseq(o7)
        catcgat
        gtagcta
        >>> c = Dseq("jatcgati")
        >>> c
        Dseq(-8)
         atcgatc
        gtagcta
        >>> c.looped()
        Dseq(o7)
        catcgat
        gtagcta
        >>> d = Dseq("ietcgazj")
        >>> d
        Dseq(-8)
        catcga
          agctag
        >>> d.looped()
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/usr/local/lib/python2.7/dist-packages/pydna/dsdna.py", line 357, in looped
            if type5 == type3 and str(sticky5) == str(rc(sticky3)):
        TypeError: DNA cannot be circularized.
        5' and 3' sticky ends not compatible!
        >>>

        """
        if self.circular:
            return copy.deepcopy(self)

        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()

        err = TypeError(
            "DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!"
        )

        if type5 != type3:
            raise err

        try:
            # Test if sticky ends are compatible
            self + self
        except TypeError:
            raise err

        new = self.cast_to_ds_left()[: len(self) - len(sticky3)]

        new.circular = True
        return new

    def five_prime_end(self) -> Tuple[str, str]:
        """Returns a 2-tuple of trings describing the structure of the 5' end of
        the DNA fragment.

        The tuple contains (type , sticky) where type is eiter "5'" or "3'".
        sticky is always in lower case and contains the sequence of the
        protruding end in 5'-3' direction.

        See examples below:


        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a = Dseq("aa", "tttg", ovhg=2)
        >>> a
        Dseq(-4)
          aa
        gttt
        >>> a.five_prime_end()
        ("3'", 'tg')
        >>> a = Dseq("caaa", "tt", ovhg=-2)
        >>> a
        Dseq(-4)
        caaa
          tt
        >>> a.five_prime_end()
        ("5'", 'ca')
        >>> a = Dseq("aa", "tt")
        >>> a
        Dseq(-2)
        aa
        tt
        >>> a.five_prime_end()
        ('blunt', '')

        See also
        --------
        pydna.dseq.Dseq.three_prime_end

        """

        # See docstring for function pydna.utils.get_parts for details
        # on what is contained in parts.
        parts = self.get_parts()

        sticky5 = parts.sticky_left5.translate(dscode_to_watson_table)

        sticky3 = parts.sticky_left3.translate(dscode_to_crick_table)[::-1]

        single_watson = parts.single_watson.translate(dscode_to_watson_table)

        single_crick = parts.single_crick.translate(dscode_to_crick_table)[::-1]

        # The walrus operator returns the value being assigned, so
        # we can test if it is empty or not.
        if sticky := single_watson:
            type_ = "single"
        elif sticky := single_crick:
            type_ = "single"
        elif sticky5 == sticky3 == "":
            type_, sticky = "blunt", ""
        elif sticky := sticky5:
            type_ = "5'"
        elif sticky := sticky3:
            type_ = "3'"

        return type_, sticky.lower()

    def three_prime_end(self) -> Tuple[str, str]:
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        >>> a = Dseq("aa", "gttt", ovhg=0)
        >>> a
        Dseq(-4)
        aa
        tttg
        >>> a.three_prime_end()
        ("5'", 'gt')
        >>> a = Dseq("aaac", "tt", ovhg=0)
        >>> a
        Dseq(-4)
        aaac
        tt
        >>> a.three_prime_end()
        ("3'", 'ac')
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', '')

        See also
        --------
        pydna.dseq.Dseq.five_prime_end

        """

        # See docstring for function pydna.utils.get_parts for details
        # on what is contained in parts.
        parts = self.get_parts()

        sticky5 = parts.sticky_right5.translate(dscode_to_crick_table)[::-1]

        sticky3 = parts.sticky_right3.translate(dscode_to_watson_table)

        single_watson = parts.single_watson.translate(dscode_to_watson_table)

        single_crick = parts.single_crick.translate(dscode_to_crick_table)[::-1]

        # The walrus operator returns the value being assigned, so
        # we can test if it is empty or not.
        if sticky := single_watson:
            type_ = "single"
        elif sticky := single_crick:
            type_ = "single"
        elif sticky5 == sticky3 == "":
            type_, sticky = "blunt", ""
        elif sticky := sticky5:
            type_ = "5'"
        elif sticky := sticky3:
            type_ = "3'"

        return type_, sticky.lower()

    def __add__(self: DseqType, other: [DseqType, str, bytes]) -> DseqType:
        """
        Adding two Dseq objects together.

        >>> ds = Dseq("a", "t", ovhg=0)
        >>> ds
        Dseq(-1)
        a
        t
        >>> ds + ds
        Dseq(-2)
        aa
        tt
        >>> "g" + ds # adding a string of left side returns a Dseq
        Dseq(-2)
        ga
        ct
        >>> ds + "c" # adding a string of right side returns a Dseq
        Dseq(-2)
        ac
        tg


        Parameters
        ----------
        other : [DseqType, str, bytes]
            Object to be added.

        Raises
        ------
        TypeError
            Preventing adding to a circular sequence.

        Returns
        -------
        DseqType
            A new Dseq object.

        """

        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        # If other evaluates to False, return a copy of self.
        if not other:
            return copy.deepcopy(self)
        # If self evaluates to False, return a copy of other.
        elif not self:
            return copy.deepcopy(other)

        # get right side end properties for self.
        self_type, self_tail = self.three_prime_end()

        try:
            other_type, other_tail = other.five_prime_end()
        except AttributeError:
            # if other does not have the expected properties
            # most likely it is a string that can be cast as
            # a Dseq.
            other_type, other_tail = "blunt", ""
            other = Dseq(other)

        err = TypeError("sticky ends not compatible!")

        # The sticky ends has to be of the same type
        # or
        # one or both of is "single" indicating a stranded molecule.
        if (self_type != other_type) and ("single" not in (self_type, other_type)):
            raise err

        # tail length has to be equal for two phosphdiester bonds to form
        if len(self_tail) != len(other_tail):
            raise err

        # Each basepair is checked against the pydna.alphabet basepair_dict
        # which contains the permitted base pairings.
        for w, c in zip(self_tail, other_tail[::-1]):
            try:
                basepair_dict[(w, c)]
            except KeyError:
                raise err

        return self.__class__(
            self.watson + other.watson, other.crick + self.crick, self.ovhg
        )

    def __mul__(self: DseqType, number: int) -> DseqType:
        if not isinstance(number, int):
            raise TypeError(
                "TypeError: can't multiply Dseq" f" by non-int of type {type(number)}"
            )
        return Dseq("").join(list(itertools.repeat(self, number)))

    def _fill_in_left(self: DseqType, nucleotides: str) -> str:
        stuffer = ""
        type, se = self.five_prime_end()
        if type == "5'":
            for n in rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.crick + stuffer, self.ovhg + len(stuffer)

    def _fill_in_right(self: DseqType, nucleotides: str) -> str:
        stuffer = ""
        type, se = self.three_prime_end()
        if type == "5'":
            for n in rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.watson + stuffer

    def fill_in(self, nucleotides: Union[None, str] = None) -> DseqType:
        """Fill in of five prime protruding end with a DNA polymerase
        that has only DNA polymerase activity (such as Exo-Klenow [#]_).
        Exo-Klenow is a modified version of the Klenow fragment of E.
        coli DNA polymerase I, which has been engineered to lack both
        3-5 proofreading and 5-3 exonuclease activities.

        and any combination of A, G, C or T. Default are all four
        nucleotides together.

        Parameters
        ----------

        nucleotides : str

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> b=Dseq("caaa", "cttt")
        >>> b
        Dseq(-5)
        caaa
         tttc
        >>> b.fill_in()
        Dseq(-5)
        caaag
        gtttc
        >>> b.fill_in("g")
        Dseq(-5)
        caaag
        gtttc
        >>> b.fill_in("tac")
        Dseq(-5)
        caaa
         tttc
        >>> c=Dseq("aaac", "tttg")
        >>> c
        Dseq(-5)
         aaac
        gttt
        >>> c.fill_in()
        Dseq(-5)
         aaac
        gttt
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.fill_in()
        Dseq(-3)
        aaa
        ttt

        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Klenow_fragment#The_exo-_Klenow_fragment

        """
        if nucleotides is None:
            nucleotides = "GATCRYWSMKHBVDN"

        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        crick, ovhg = self._fill_in_left(nucleotides)
        watson = self._fill_in_right(nucleotides)
        return Dseq(watson, crick, ovhg)

    klenow = fill_in  # alias

    def nibble_to_blunt(self) -> DseqType:
        """
        Simulates treatment a nuclease with both 5'-3' and 3'-5' single
        strand specific exonuclease activity (such as mung bean nuclease [#]_)

        Mung bean nuclease is a nuclease enzyme derived from mung bean sprouts
        that preferentially degrades single-stranded DNA and RNA into
        5'-phosphate- and 3'-hydroxyl-containing nucleotides.

        Treatment results in blunt DNA, regardless of wheter the protruding end
        is 5' or 3'.

        ::

             ggatcc    ->     gatcc
              ctaggg          ctagg

              ggatcc   ->     ggatc
             tcctag           cctag

         >>> from pydna.dseq import Dseq
         >>> b=Dseq("caaa", "cttt")
         >>> b
         Dseq(-5)
         caaa
          tttc
         >>> b.mung()
         Dseq(-3)
         aaa
         ttt
         >>> c=Dseq("aaac", "tttg")
         >>> c
         Dseq(-5)
          aaac
         gttt
         >>> c.mung()
         Dseq(-3)
         aaa
         ttt



        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Mung_bean_nuclease


        """
        parts = self.get_parts()
        return self.__class__(parts.middle)

    mung = nibble_to_blunt

    def T4(self, nucleotides=None) -> DseqType:
        """
        Fill in 5' protruding ends and nibble 3' protruding ends.

        This is done using a DNA polymerase providing 3'-5' nuclease activity
        such as T4 DNA polymerase. This can be done in presence of any
        combination of the four nucleotides A, G, C or T.

        T4 DNA polymerase is widely used to “polish” DNA ends because of its
        strong 3-5 exonuclease activity in the absence of dNTPs, it chews
        back 3′ overhangs to create blunt ends; in the presence of limiting
        dNTPs, it can fill in 5′ overhangs; and by carefully controlling
        reaction time, temperature, and nucleotide supply, you can generate
        defined recessed or blunt termini.

        Tuning the nucleotide set can facilitate engineering of partial
        sticky ends. Default are all four nucleotides together.

        ::

                  aaagatc-3        aaa      3' ends are always removed.
                  |||       --->   |||      A and T needed or the molecule will
            3-ctagttt              ttt      degrade completely.



            5-gatcaaa              gatcaaaGATC      5' ends are filled in the
                  |||       --->   |||||||||||      presence of GATC
                  tttctag-5        CTAGtttctag



            5-gatcaaa              gatcaaaGAT       5' ends are partially filled in the
                  |||       --->    |||||||||       presence of GAT to produce a 1 nt
                  tttctag-5         TAGtttctag      5' overhang



            5-gatcaaa              gatcaaaGA       5' ends are partially filled in the
                  |||       --->     |||||||       presence of GA to produce a 2 nt
                  tttctag-5          AGtttctag     5' overhang



            5-gatcaaa              gatcaaaG        5' ends are partially filled in the
                  |||       --->      |||||        presence of G to produce a 3 nt
                  tttctag-5           Gtttctag     5' overhang



        Parameters
        ----------
        nucleotides : str


        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a = Dseq.from_representation(
        ... '''
        ... gatcaaa
        ...     tttctag
        ... ''')
        >>> a
        Dseq(-11)
        gatcaaa
            tttctag
        >>> a.T4()
        Dseq(-11)
        gatcaaagatc
        ctagtttctag
        >>> a.T4("GAT")
        Dseq(-11)
        gatcaaagat
         tagtttctag
        >>> a.T4("GA")
        Dseq(-11)
        gatcaaaga
          agtttctag
        >>> a.T4("G")
        Dseq(-11)
        gatcaaag
           gtttctag
        """

        if not nucleotides:
            nucleotides = "GATCRYWSMKHBVDN"
        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        type, se = self.five_prime_end()
        if type == "5'":
            crick, ovhg = self._fill_in_left(nucleotides)
        else:
            if type == "3'":
                ovhg = 0
                crick = self.crick[: -len(se)]
            else:
                ovhg = 0
                crick = self.crick
        x = len(crick) - 1
        while x >= 0:
            if crick[x] in nucleotides:
                break
            x -= 1
        ovhg = x - len(crick) + 1 + ovhg
        crick = crick[: x + 1]
        if not crick:
            ovhg = 0
        watson = self.watson
        type, se = self.three_prime_end()
        if type == "5'":
            watson = self._fill_in_right(nucleotides)
        else:
            if type == "3'":
                watson = self.watson[: -len(se)]
        x = len(watson) - 1
        while x >= 0:
            if watson[x] in nucleotides:
                break
            x -= 1
        watson = watson[: x + 1]
        return Dseq(watson, crick, ovhg)

    t4 = T4  # alias for the T4 method.

    def nibble_five_prime_left(self: DseqType, n: int = 1) -> DseqType:
        """
        5' => 3'  resection at the left side (start) of the molecule.

        The argument n indicate the number of nucleotides that are to be
        removed. The outcome of this depend on the structure of the molecule.
        See the two examples below:

        The figure below indicates a recess of length two from a blunt DNA
        fragment. The resulting DNA fragment has a 3' protruding single strand.

        ::

            gatc           tc
            ||||   -->     ||
            ctag         ctag


        The figure below indicates a recess of length two from a DNA fragment
        with a 5' sticky end resulting in a blunt sequence.

        ::

          ttgatc         gatc
            ||||   -->   ||||
            ctag         ctag


        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("gatc")
        >>> ds
        Dseq(-4)
        gatc
        ctag
        >>> ds.nibble_five_prime_left(2)
        Dseq(-4)
          tc
        ctag
        >>> ds.nibble_five_prime_left(3)
        Dseq(-4)
           c
        ctag
        >>> ds.nibble_five_prime_left(4)
        Dseq(-4)
        <BLANKLINE>
        ctag
        >>> ds = Dseq.from_representation(
        ... '''
        ... GGgatc
        ...   ctag
        ... ''')
        >>> ds
        Dseq(-6)
        GGgatc
          ctag
        >>> ds.nibble_five_prime_left(2)
        Dseq(-4)
        gatc
        ctag

        Parameters
        ----------
        n : int, optional
            The default is 1. This is the number of nucleotides removed.

        Returns
        -------
        DseqType
            DESCRIPTION.

        """
        n += max(0, self.ovhg or 0)
        return Dseq(
            self._data[:n]
            .translate(dscode_to_crick_table)
            .translate(complement_table_for_dscode)
            .translate(dscode_to_crick_tail_table)
            .lstrip()
            + self._data[n:]
        )

    def nibble_five_prime_right(self: DseqType, n: int = 1) -> DseqType:
        """
        5' => 3'  resection at the right side (end) of the molecule.

        The argument n indicate the number of nucleotides that are to be
        removed. The outcome of this depend on the structure of the molecule.
        See the two examples below:

        The figure below indicates a recess of length two from a blunt DNA
        fragment. The resulting DNA fragment has a 3' protruding single strand.

        ::

            gatc         gatc
            ||||   -->   ||
            ctag         ct

        The figure below indicates a recess of length two from a DNA fragment
        with a 5' sticky end resulting in a blunt sequence.

        ::

            gatc         gatc
            ||||   -->   ||||
            ctagtt       ctag


        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("gatc")
        >>> ds
        Dseq(-4)
        gatc
        ctag
        >>> ds.nibble_five_prime_right(2)
        Dseq(-4)
        gatc
        ct
        >>> ds.nibble_five_prime_right(3)
        Dseq(-4)
        gatc
        c
        >>> ds.nibble_five_prime_right(4)
        Dseq(-4)
        gatc
        <BLANKLINE>
        >>> ds = Dseq.from_representation(
        ... '''
        ... gatc
        ... ctagGG
        ... ''')
        >>> ds.nibble_five_prime_right(2)
        Dseq(-4)
        gatc
        ctag
        """
        n = len(self) - n
        ovhg = len(self) if self.right_ovhg is None else self.right_ovhg
        n -= max(0, ovhg)
        return Dseq(
            self._data[:n]
            + self._data[n:]
            .translate(dscode_to_watson_table)
            .translate(dscode_to_watson_tail_table)
            .lstrip()
        )

    exo1_front = nibble_five_prime_left  # TODO: consider using the new names
    exo1_end = nibble_five_prime_right  # TODO: consider using the new names

    def nibble_three_prime_left(self: DseqType, n=1) -> DseqType:
        """
        3' => 5' resection at the left side (beginning) of the molecule.

        The argument n indicate the number of nucleotides that are to be
        removed. The outcome of this depend on the structure of the molecule.
        See the two examples below:

        The figure below indicates a recess of length two from a blunt DNA
        fragment. The resulting DNA fragment has a 5' protruding single strand.

        ::

            gatc         gatc
            ||||   -->     ||
            ctag           ag

        The figure below indicates a recess of length two from a DNA fragment
        with a 3' sticky end resulting in a blunt sequence.

        ::

            gatc         gatc
            ||||   -->   ||||
          ttctag         ctag


        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("gatc")
        >>> ds
        Dseq(-4)
        gatc
        ctag
        >>> ds.nibble_three_prime_left(2)
        Dseq(-4)
        gatc
          ag
        >>> ds.nibble_three_prime_left(3)
        Dseq(-4)
        gatc
           g
        >>> ds.nibble_three_prime_left(4)
        Dseq(-4)
        gatc
        <BLANKLINE>
        >>> ds = Dseq.from_representation(
        ... '''
        ...   gatc
        ... CCctag
        ... ''')
        >>> ds
        Dseq(-6)
          gatc
        CCctag
        >>> ds.nibble_three_prime_left(2)
        Dseq(-4)
        gatc
        ctag
        """
        ovhg = len(self) if self.ovhg is None else self.ovhg
        n -= min(0, ovhg)
        return Dseq(
            self._data[:n]
            .translate(dscode_to_watson_table)
            .translate(dscode_to_watson_tail_table)
            .lstrip()
            + self._data[n:]
        )

    def nibble_three_prime_right(self: DseqType, n=1) -> DseqType:
        """
        3' => 5' resection at the right side (end) of the molecule.

        The argument n indicate the number of nucleotides that are to be
        removed. The outcome of this depend on the structure of the molecule.
        See the two examples below:

        The figure below indicates a recess of length two from a blunt DNA
        fragment. The resulting DNA fragment has a 5' protruding single strand.

        ::

            gatc         ga
            ||||   -->   ||
            ctag         ctag

        The figure below indicates a recess of length two from a DNA fragment
        with a 3' sticky end resulting in a blunt sequence.

        ::

            gatctt       gatc
            ||||   -->   ||||
            ctag         ctag


        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("gatc")
        >>> ds
        Dseq(-4)
        gatc
        ctag
        >>> ds.nibble_three_prime_right(2)
        Dseq(-4)
        ga
        ctag
        >>> ds.nibble_three_prime_right(3)
        Dseq(-4)
        g
        ctag
        >>> ds.nibble_three_prime_right(4)
        Dseq(-4)
        <BLANKLINE>
        ctag
        >>> ds = Dseq.from_representation(
        ... '''
        ... gatcCC
        ... ctag
        ... ''')
        >>> ds.nibble_three_prime_right(2)
        Dseq(-4)
        gatc
        ctag
        """
        n = len(self) - n
        ovhg = len(self) if self.right_ovhg is None else self.right_ovhg
        n += min(0, ovhg)
        return Dseq(
            self._data[:n]
            + self._data[n:]
            .translate(dscode_to_crick_table)
            .translate(complement_table_for_dscode)
            .translate(dscode_to_crick_tail_table)
            .lstrip()
        )

    def no_cutters(
        self, batch: Union[RestrictionBatch, None] = None
    ) -> RestrictionBatch:
        """Enzymes in a RestrictionBatch not cutting sequence."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if not sitelist}
        return RestrictionBatch(ncut)

    def unique_cutters(
        self, batch: Union[RestrictionBatch, None] = None
    ) -> RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence once."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=1, batch=batch)

    once_cutters = unique_cutters  # alias for unique_cutters

    def twice_cutters(
        self, batch: Union[RestrictionBatch, None] = None
    ) -> RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence twice."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=2, batch=batch)

    def n_cutters(
        self, n=3, batch: Union[RestrictionBatch, None] = None
    ) -> RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting n times."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if len(sitelist) == n}
        return RestrictionBatch(ncut)

    def cutters(self, batch: Union[RestrictionBatch, None] = None) -> RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence at least once."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if sitelist}
        return RestrictionBatch(ncut)

    def seguid(self) -> str:
        """SEGUID checksum for the sequence."""
        if self.circular:
            cs = cdseguid(
                self.watson.upper(), self.crick.upper(), alphabet="{DNA-extended}"
            )
        else:
            """docstring."""
            w = f"{self.ovhg * '-'}{self.watson}{'-' * (-self.ovhg + len(self.crick) - len(self.watson))}".upper()
            c = f"{'-' * (self.ovhg + len(self.watson) - len(self.crick))}{self.crick}{-self.ovhg * '-'}".upper()
            cs = ldseguid(w, c, alphabet="{DNA-extended},AU")
        return cs

    def isblunt(self) -> bool:
        """isblunt.

        Return True if Dseq is linear and blunt and
        false if staggered or circular.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("gat")
        >>> a
        Dseq(-3)
        gat
        cta
        >>> a.isblunt()
        True
        >>> a=Dseq("gat", "atcg")
        >>> a
        Dseq(-4)
         gat
        gcta
        >>> a.isblunt()
        False
        >>> a=Dseq("gat", "gatc")
        >>> a
        Dseq(-4)
        gat
        ctag
        >>> a.isblunt()
        False
        >>> a=Dseq("gat", circular=True)
        >>> a
        Dseq(o3)
        gat
        cta
        >>> a.isblunt()
        False
        """
        parts = self.get_parts()

        return not any(
            (
                parts.sticky_right5,
                parts.sticky_right3,
                parts.sticky_left3,
                parts.sticky_left5,
                self.circular,
            )
        )

    def terminal_transferase(self, nucleotides: str = "a") -> DseqType:
        """
        Terminal deoxynucleotidyl transferase (TdT) is a template-independent
        DNA polymerase that adds nucleotides to the 3′-OH ends of DNA, typically
        single-stranded or recessed 3′ ends. In cloning, it’s classically used
        to create homopolymer tails (e.g. poly-dG on a vector and poly-dC on an insert)
        so that fragments can anneal via complementary overhangs (“tailing” cloning).

        This activity ia also present in some DNA polymerases, such as Taq polymerase.
        This property is used in the populat T/A cloning protocol ([#]_).

        ::

            gct          gcta
            |||   -->    |||
            cga         acga



        >>> from pydna.dseq import Dseq
        >>> a = Dseq("aa")
        >>> a = Dseq("gct")
        >>> a
        Dseq(-3)
        gct
        cga
        >>> a.terminal_transferase()
        Dseq(-5)
         gcta
        acga
        >>> a.terminal_transferase("G")
        Dseq(-5)
         gctG
        Gcga

        Parameters
        ----------
        nucleotides : str, optional
            The default is "a".

        Returns
        -------
        DseqType
            DESCRIPTION.

        References
        ----------
        .. [#] https://en.wikipedia.org/wiki/TA_cloning

        """
        ovhg = self.ovhg
        if self.ovhg >= 0:
            ovhg += len(nucleotides)
        return Dseq(self.watson + nucleotides, self.crick + nucleotides, ovhg)

    def user(self) -> DseqType:
        """
        USER Enzyme treatment.

        USER Enzyme is a mixture of Uracil DNA glycosylase (UDG) and the
        DNA glycosylase-lyase Endonuclease VIII.

        UDG catalyses the excision of an uracil base, forming an abasic
        or apyrimidinic site (AP site). Endonuclease VIII removes the AP
        site creating a DNA gap.

        ::

            tagaagtaggUat          tagaagtagg at
            |||||||||||||  --->    |||||||||| ||
            atcUtcatccata          atc tcatccata



        >>> a = Dseq("tagaagtaggUat", "atcUtcatccata"[::-1], 0)
        >>> a
        Dseq(-13)
        tagaagtaggUat
        atcutcatccAta
        >>> a.user()
        Dseq(-13)
        tagaagtagg at
        atc tcatccAta


        Returns
        -------
        DseqType
            DNA fragment with uracile bases removed.

        """

        return Dseq(self._data.translate(bytes.maketrans(b"UuOo", b"ZzEe")))

    def cut(self: DseqType, *enzymes: EnzymesType) -> Tuple[DseqType, ...]:
        """Returns a list of linear Dseq fragments produced in the digestion.
        If there are no cuts, an empty list is returned.

        Parameters
        ----------

        enzymes : enzyme object or iterable of such objects
            A Bio.Restriction.XXX restriction objects or iterable.

        Returns
        -------
        frags : list
            list of Dseq objects formed by the digestion


        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> seq=Dseq("ggatccnnngaattc")
        >>> seq
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>> from Bio.Restriction import BamHI,EcoRI
        >>> type(seq.cut(BamHI))
        <class 'tuple'>
        >>> for frag in seq.cut(BamHI): print(repr(frag))
        Dseq(-5)
        g
        cctag
        Dseq(-14)
        gatccnnngaattc
            gnnncttaag
        >>> seq.cut(EcoRI, BamHI) ==  seq.cut(BamHI, EcoRI)
        True
        >>> a,b,c = seq.cut(EcoRI, BamHI)
        >>> a+b+c
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>>

        """

        cutsites = self.get_cutsites(*enzymes)
        cutsite_pairs = self.get_cutsite_pairs(cutsites)
        return tuple(self.apply_cut(*cs) for cs in cutsite_pairs)

    def cutsite_is_valid(self, cutsite: CutSiteType) -> bool:
        """
        Check is a cutsite is valid.

        A cutsite is a nested 2-tuple with this form:

        ((cut_watson, ovhg), enz), for example ((396, -4), EcoRI)

        The cut_watson (positive integer) is the cut position of the sequence as for example
        returned by the Bio.Restriction module.

        The ovhg (overhang, positive or negative integer or 0) has the same meaning as
        for restriction enzymes in the Bio.Restriction module and for
        pydna.dseq.Dseq objects (see docstring for this module and example below)

        Enzyme can be None.

        ::

            Enzyme overhang

            EcoRI  -4     --GAATTC--        --G       AATTC--
                            ||||||     -->    |           |
                          --CTTAAG--        --CTTAA       G--

            KpnI    4     --GGTACC--        --GGTAC       C--
                            ||||||     -->    |           |
                          --CCATGG--        --C       CATGG--

            SmaI    0     --CCCGGG--        --CCC       GGG--
                            ||||||     -->    |||       |||
                          --GGGCCC--        --GGG       CCC--


        >>> from Bio.Restriction import EcoRI, KpnI, SmaI
        >>> EcoRI.ovhg
        -4
        >>> KpnI.ovhg
        4
        >>> SmaI.ovhg
        0

        Returns False if:

        - Cut positions fall outside the sequence (could be moved to Biopython)
        TODO: example

        - Overhang is not double stranded
        TODO: example

        - Recognition site is not double stranded or is outside the sequence
        TODO: example

        - For enzymes that cut twice, it checks that at least one possibility is valid
        TODO: example



        Parameters
        ----------
        cutsite : CutSiteType
            DESCRIPTION.

        Returns
        -------
        bool
            True if cutsite can cut the DNA fragment.

        """

        assert cutsite is not None, "cutsite is None"

        enz = cutsite[1]
        watson, crick, ovhg = self.get_cut_parameters(cutsite, True)

        # The overhang is double stranded
        overhang_dseq = self[watson:crick] if ovhg < 0 else self[crick:watson]
        if overhang_dseq.ovhg != 0 or overhang_dseq.watson_ovhg != 0:
            return False

        # The recognition site is double stranded and within the sequence
        start_of_recognition_site = watson - enz.fst5
        if start_of_recognition_site < 0:
            start_of_recognition_site += len(self)
        end_of_recognition_site = start_of_recognition_site + enz.size
        if self.circular:
            end_of_recognition_site %= len(self)
        recognition_site = self[start_of_recognition_site:end_of_recognition_site]
        if (
            len(recognition_site) == 0
            or recognition_site.ovhg != 0
            or recognition_site.watson_ovhg != 0
        ):
            if enz is None or enz.scd5 is None:
                return False
            else:
                # For enzymes that cut twice, this might be referring to the second one
                start_of_recognition_site = watson - enz.scd5
                if start_of_recognition_site < 0:
                    start_of_recognition_site += len(self)
                end_of_recognition_site = start_of_recognition_site + enz.size
                if self.circular:
                    end_of_recognition_site %= len(self)
                recognition_site = self[
                    start_of_recognition_site:end_of_recognition_site
                ]

                if (
                    len(recognition_site) == 0
                    or recognition_site.ovhg != 0
                    or recognition_site.watson_ovhg != 0
                ):
                    return False

        return True

    def get_cutsites(self: DseqType, *enzymes: EnzymesType) -> List[CutSiteType]:
        """Returns a list of cutsites, represented represented as `((cut_watson, ovhg), enz)`:

        - `cut_watson` is a positive integer contained in `[0,len(seq))`, where `seq` is the sequence
          that will be cut. It represents the position of the cut on the watson strand, using the full
          sequence as a reference. By "full sequence" I mean the one you would get from `str(Dseq)`.

        - `ovhg` is the overhang left after the cut. It has the same meaning as `ovhg` in
          the `Bio.Restriction` enzyme objects, or pydna's `Dseq` property.

        - `enz` is the enzyme object. It's not necessary to perform the cut, but can be
           used to keep track of which enzyme was used.

        Cuts are only returned if the recognition site and overhang are on the double-strand
        part of the sequence.

        Parameters
        ----------

        enzymes : Union[RestrictionBatch,list[_AbstractCut]]

        Returns
        -------
        list[tuple[tuple[int,int], _AbstractCut]]

        Examples
        --------

        >>> from Bio.Restriction import EcoRI
        >>> from pydna.dseq import Dseq
        >>> seq = Dseq('AAGAATTCAAGAATTC')
        >>> seq.get_cutsites(EcoRI)
        [((3, -4), EcoRI), ((11, -4), EcoRI)]

        `cut_watson` is defined with respect to the "full sequence", not the
        watson strand:

        >>> dseq = Dseq.from_full_sequence_and_overhangs('aaGAATTCaa', 1, 0)
        >>> dseq
        Dseq(-10)
         aGAATTCaa
        ttCTTAAGtt
        >>> dseq.get_cutsites([EcoRI])
        [((3, -4), EcoRI)]

        Cuts are only returned if the recognition site and overhang are on the double-strand
        part of the sequence.

        >>> Dseq('GAATTC').get_cutsites([EcoRI])
        [((1, -4), EcoRI)]
        >>> Dseq.from_full_sequence_and_overhangs('GAATTC', -1, 0).get_cutsites([EcoRI])
        []

        """

        if len(enzymes) == 1 and isinstance(enzymes[0], RestrictionBatch):
            # argument is probably a RestrictionBatch
            enzymes = [e for e in enzymes[0]]

        enzymes = list(dict.fromkeys(flatten(enzymes)))  # remove duplicate enzymes
        out = list()
        for e in enzymes:
            # Positions of the cut on the watson strand. They are 1-based, so we subtract
            # 1 to get 0-based positions
            cuts_watson = [c - 1 for c in e.search(self, linear=(not self.circular))]

            out += [((w, e.ovhg), e) for w in cuts_watson]

        return sorted([cutsite for cutsite in out if self.cutsite_is_valid(cutsite)])

    def left_end_position(self) -> Tuple[int, int]:
        """
        The index in the full sequence of the watson and crick start positions.

        full sequence (str(self)) for all three cases is AAA
        ::

            AAA              AA               AAT
             TT             TTT               TTT
            Returns (0, 1)  Returns (1, 0)    Returns (0, 0)


        """
        if self.ovhg > 0:
            return self.ovhg, 0
        return 0, -self.ovhg

    def right_end_position(self) -> Tuple[int, int]:
        """The index in the full sequence of the watson and crick end positions.

        full sequence (str(self)) for all three cases is AAA

        ```
        AAA               AA                   AAA
        TT                TTT                  TTT
        Returns (3, 2)    Returns (2, 3)       Returns (3, 3)
        ```

        """
        if self.watson_ovhg < 0:
            return len(self) + self.watson_ovhg, len(self)
        return len(self), len(self) - self.watson_ovhg

    def get_ss_meltsites(self: DseqType, length: int) -> tuple[int, int]:
        """
        Single stranded DNA melt sites

        Two lists of 2-tuples of integers are returned. Each tuple
        (`((from, to))`) contains the start and end positions of a single
        stranded region, shorter or equal to `length`.

        In the example below, the middle 2 nt part is released from the
        molecule.

        ::


            tagaa ta gtatg
            ||||| || |||||  -->   [(6,8)], []
            atcttcatccatac

            tagaagtaggtatg
            ||||| || |||||  -->   [], [(6,8)]
            atctt at catac




        The output of this method is used in the `melt_ss_dna` method in order
        to determine the start and end positions of single stranded regions.

        See get_ds_meltsites for melting ds sequences.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("tagaaqtaqgtatg")
        >>> ds
        Dseq(-14)
        tagaa ta gtatg
        atcttcatccatac
        >>> cutsites = ds.get_ss_meltsites(2)
        >>> cutsites
        ([(6, 8)], [])
        >>> ds[6:8]
        Dseq(-2)
        ta
        at
        >>> ds = Dseq("tagaaptapgtatg")
        >>> ds
        Dseq(-14)
        tagaagtaggtatg
        atctt at catac
        >>> cutsites = ds.get_ss_meltsites(2)
        >>> cutsites
        ([], [(6, 8)])
        """

        regex = regex_ss_melt_factory(length)

        if self.circular:
            spacer = length
            cutfrom = self._data[-length:] + self._data + self._data[:length]
        else:
            spacer = 0
            cutfrom = self._data

        watson_cuts = []
        crick_cuts = []

        for m in regex.finditer(cutfrom):

            if m.lastgroup == "watson":
                cut1 = m.start() + spacer
                cut2 = m.end() + spacer
                watson_cuts.append((cut1, cut2))
            else:
                assert m.lastgroup == "crick"
                cut1 = m.start() + spacer
                cut2 = m.end() + spacer
                crick_cuts.append((cut1, cut2))

        return watson_cuts, crick_cuts

    def get_ds_meltsites(self: DseqType, length: int) -> List[CutSiteType]:
        """
        Double stranded DNA melt sites

        DNA molecules can fall apart by melting if they have internal single
        stranded regions. In the example below, the molecule has two gaps
        on opposite sides, two nucleotides apart, which means that it hangs
        together by two basepairs.

        This molecule can melt into two separate 8 bp double stranded
        molecules, each with 3 nt 3' overhangs a depicted below.

        ::

            tagaagta gtatg        tagaagta          gtatg
            ||||| || |||||  -->   |||||             |||||
            atctt atccatac        atctt          atccatac


        A list of 2-tuples is returned. Each tuple (`((cut_watson, ovhg), None)`)
        contains cut position and the overhang value in the same format as
        returned by the get_cutsites method for restriction enzymes.

        Note that this function deals with melting that results in two double
        stranded DNA molecules.

        See get_ss_meltsites for melting of single stranded regions from
        molecules.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("tagaaptaqgtatg")
        >>> ds
        Dseq(-14)
        tagaagta gtatg
        atctt atccatac
        >>> cutsite = ds.get_ds_meltsites(2)
        >>> cutsite
        [((8, 2), None)]

        """

        if length < 1:
            return tuple()

        regex = regex_ds_melt_factory(length)

        if self.circular:
            spacer = length
            cutfrom = self._data[-length:] + self._data + self._data[:length]
        else:
            spacer = 0
            cutfrom = self._data

        cuts = []

        for m in regex.finditer(cutfrom):

            if m.lastgroup == "watson":
                cut = (m.end() - spacer, m.end() - m.start()), None
            else:
                assert m.lastgroup == "crick"
                cut = (m.start() - spacer, m.start() - m.end()), None

            cuts.append(cut)

        return cuts

    def cast_to_ds_right(self):
        """
        NNNN               NNNNGATC
        ||||       -->     ||||||||
        NNNNCTAG           NNNNCTAG


        NNNNGATC           NNNNGATC
        ||||       -->     ||||||||
        NNNN               NNNNCTAG
        """

        p = self.get_parts()

        ds_stuffer = (p.sticky_right5 or p.sticky_right3).translate(
            dscode_to_full_sequence_table
        )

        result = (p.sticky_left5 or p.sticky_left3) + p.middle + ds_stuffer

        return self.__class__(result, circular=False)

    def cast_to_ds(self):
        """Sequencially calls cast_to_ds_left and cast_to_ds_right."""
        return self.cast_to_ds_left().cast_to_ds_right()

    def cast_to_ds_left(self):
        """
        GATCNNNN           GATCNNNN
            ||||   -->     ||||||||
            NNNN           CTAGNNNN

            NNNN           GATCNNNN
            ||||   -->     ||||||||
        CTAGNNNN           CTAGNNNN
        """

        p = self.get_parts()

        ds_stuffer = (p.sticky_left5 or p.sticky_left3).translate(
            dscode_to_full_sequence_table
        )

        result = ds_stuffer + p.middle + (p.sticky_right5 or p.sticky_right3)

        return self.__class__(result, circular=False)

    def get_cut_parameters(
        self, cut: Union[CutSiteType, None], is_left: bool
    ) -> Tuple[int, int, int]:
        """For a given cut expressed as ((cut_watson, ovhg), enz), returns
        a tuple (cut_watson, cut_crick, ovhg).

        - cut_watson: see get_cutsites docs
        - cut_crick: equivalent of cut_watson in the crick strand
        - ovhg: see get_cutsites docs

        The cut can be None if it represents the left or right end of the sequence.
        Then it will return the position of the watson and crick ends with respect
        to the "full sequence". The `is_left` parameter is only used in this case.

        """
        if cut is not None:
            watson, ovhg = cut[0]
            crick = watson - ovhg
            if self.circular:
                crick %= len(self)
            return watson, crick, ovhg

        assert not self.circular, "Circular sequences should not have None cuts"

        if is_left:
            return *self.left_end_position(), self.ovhg
        # In the right end, the overhang does not matter
        return *self.right_end_position(), self.watson_ovhg

    def melt(self, length):
        """
        TBD

        Parameters
        ----------
        length : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if not length or length < 1:
            return tuple()

        # First we need to get rid of single stranded sequences
        new, strands = self.melt_ss_dna(length)

        cutsites = new.get_ds_meltsites(length)

        cutsite_pairs = self.get_cutsite_pairs(cutsites)

        result = tuple(new.apply_cut(*cutsite_pair) for cutsite_pair in cutsite_pairs)

        result = tuple([new]) if strands and not result else result

        return tuple(strands) + tuple(result)

    def melt_ss_dna(self, length) -> tuple["Dseq", list["Dseq"]]:
        """
        Melt to separate single stranded DNA

        Single stranded DNA molecules shorter or equal to `length` shed from
        a double stranded DNA molecule without affecting the length of the
        remaining molecule.

        In the examples below, the middle 2 nt part is released from the
        molecule.

        ::

            tagaa ta gtatg        tagaa    gtatg          ta
            ||||| || |||||  -->   |||||    |||||     +    ||
            atcttcatccatac        atcttcatccatac

            tagaagtaggtatg        tagaagtaggtatg
            ||||| || |||||  -->   |||||    |||||     +    ||
            atctt at catac        atctt    catac          at


        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("tagaaqtaqgtatg")
        >>> ds
        Dseq(-14)
        tagaa ta gtatg
        atcttcatccatac
        >>> new, strands  = ds.melt_ss_dna(2)
        >>> new
        Dseq(-14)
        tagaa    gtatg
        atcttcatccatac
        >>> strands[0]
        Dseq(-2)
        ta
        <BLANKLINE>
        >>> ds = Dseq("tagaaptapgtatg")
        >>> ds
        Dseq(-14)
        tagaagtaggtatg
        atctt at catac
        >>> new, strands = ds.melt_ss_dna(2)
        >>> new
        Dseq(-14)
        tagaagtaggtatg
        atctt    catac
        >>> strands[0]
        Dseq(-2)
        <BLANKLINE>
        at
        """

        watsonnicks, cricknicks = self.get_ss_meltsites(length)

        new, strands = self.shed_ss_dna(watsonnicks, cricknicks)

        return new, strands

    def shed_ss_dna(
        self,
        watson_cutpairs: list[tuple[int, int]] = None,
        crick_cutpairs: list[tuple[int, int]] = None,
    ):
        """
        Separate parts of one of the DNA strands

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("tagaagtaggtatg")
        >>> ds
        Dseq(-14)
        tagaagtaggtatg
        atcttcatccatac
        >>> new, strands = ds.shed_ss_dna([(6, 8)],[])
        >>> new
        Dseq(-14)
        tagaag  ggtatg
        atcttcatccatac
        >>> strands[0]
        Dseq(-2)
        ta
        <BLANKLINE>
        >>> new, strands = ds.shed_ss_dna([],[(6, 8)])
        >>> new
        Dseq(-14)
        tagaagtaggtatg
        atcttc  ccatac
        >>> strands[0]
        Dseq(-2)
        <BLANKLINE>
        at
        >>> ds = Dseq("tagaagtaggtatg")
        >>> new, (strand1, strand2) = ds.shed_ss_dna([(6, 8), (9, 11)],[])
        >>> new
        Dseq(-14)
        tagaag  g  atg
        atcttcatccatac
        >>> strand1
        Dseq(-2)
        ta
        <BLANKLINE>
        >>> strand2
        Dseq(-2)
        gt
        <BLANKLINE>
        """

        watson_cutpairs = watson_cutpairs or list()
        crick_cutpairs = crick_cutpairs or list()
        strands = []

        new = bytearray(self._data)

        for x, y in watson_cutpairs:
            stuffer = new[x:y]
            ss = Dseq.quick(new[x:y].translate(dscode_to_watson_tail_table))
            new[x:y] = stuffer.translate(dscode_to_crick_tail_table)
            strands.append(ss)

        for x, y in crick_cutpairs:
            stuffer = new[x:y]
            ss = Dseq.quick(stuffer.translate(dscode_to_crick_tail_table))
            new[x:y] = stuffer.translate(dscode_to_watson_tail_table)
            strands.append(ss)

        return Dseq.quick(new), strands

    def apply_cut(self, left_cut: CutSiteType, right_cut: CutSiteType) -> "Dseq":
        """Extracts a subfragment of the sequence between two cuts.

        For more detail see the documentation of get_cutsite_pairs.

        Parameters
        ----------
        left_cut : Union[tuple[tuple[int,int], _AbstractCut], None]
        right_cut: Union[tuple[tuple[int,int], _AbstractCut], None]

        Returns
        -------
        Dseq

        Examples
        --------
        >>> from Bio.Restriction import EcoRI
        >>> from pydna.dseq import Dseq
        >>> dseq = Dseq('aaGAATTCaaGAATTCaa')
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((3, -4), EcoRI), ((11, -4), EcoRI)]
        >>> p1, p2, p3 = dseq.get_cutsite_pairs(cutsites)
        >>> p1
        (None, ((3, -4), EcoRI))
        >>> dseq.apply_cut(*p1)
        Dseq(-7)
        aaG
        ttCTTAA
        >>> p2
        (((3, -4), EcoRI), ((11, -4), EcoRI))
        >>> dseq.apply_cut(*p2)
        Dseq(-12)
        AATTCaaG
            GttCTTAA
        >>> p3
        (((11, -4), EcoRI), None)
        >>> dseq.apply_cut(*p3)
        Dseq(-7)
        AATTCaa
            Gtt

        >>> dseq = Dseq('TTCaaGAA', circular=True)
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((6, -4), EcoRI)]
        >>> pair = dseq.get_cutsite_pairs(cutsites)[0]
        >>> pair
        (((6, -4), EcoRI), ((6, -4), EcoRI))
        >>> dseq.apply_cut(*pair)
        Dseq(-12)
        AATTCaaG
            GttCTTAA

        """
        if cuts_overlap(left_cut, right_cut, len(self)):
            raise ValueError("Cuts by {} {} overlap.".format(left_cut[1], right_cut[1]))

        left_watson, left_crick, ovhg_left = self.get_cut_parameters(left_cut, True)
        right_watson, right_crick, _ = self.get_cut_parameters(right_cut, False)
        return Dseq(
            self[left_watson:right_watson]._data.translate(dscode_to_watson_table),
            self[left_crick:right_crick]
            .reverse_complement()
            ._data.translate(dscode_to_watson_table),
            ovhg=ovhg_left,
        )

    def get_cutsite_pairs(
        self, cutsites: List[CutSiteType]
    ) -> List[Tuple[Union[None, CutSiteType], Union[None, CutSiteType]]]:
        """Returns pairs of cutsites that render the edges of the resulting fragments.

        A fragment produced by restriction is represented by a tuple of length 2 that
        may contain cutsites or `None`:

            - Two cutsites: represents the extraction of a fragment between those two
              cutsites, in that orientation. To represent the opening of a circular
              molecule with a single cutsite, we put the same cutsite twice.
            - `None`, cutsite: represents the extraction of a fragment between the left
              edge of linear sequence and the cutsite.
            - cutsite, `None`: represents the extraction of a fragment between the cutsite
              and the right edge of a linear sequence.

        Parameters
        ----------
        cutsites : list[tuple[tuple[int,int], _AbstractCut]]

        Returns
        -------
        list[tuple[tuple[tuple[int,int], _AbstractCut]|None],tuple[tuple[int,int], _AbstractCut]|None]

        Examples
        --------

        >>> from Bio.Restriction import EcoRI
        >>> from pydna.dseq import Dseq
        >>> dseq = Dseq('aaGAATTCaaGAATTCaa')
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((3, -4), EcoRI), ((11, -4), EcoRI)]
        >>> dseq.get_cutsite_pairs(cutsites)
        [(None, ((3, -4), EcoRI)), (((3, -4), EcoRI), ((11, -4), EcoRI)), (((11, -4), EcoRI), None)]

        >>> dseq = Dseq('TTCaaGAA', circular=True)
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((6, -4), EcoRI)]
        >>> dseq.get_cutsite_pairs(cutsites)
        [(((6, -4), EcoRI), ((6, -4), EcoRI))]
        """
        if len(cutsites) == 0:
            return []
        if not self.circular:
            cutsites = [None, *cutsites, None]
        else:
            # Add the first cutsite at the end, for circular cuts
            cutsites.append(cutsites[0])

        return list(zip(cutsites, cutsites[1:]))

    def get_parts(self):
        """
        Returns a DseqParts instance containing the parts (strings) of a dsDNA
        sequence. DseqParts instance field names:

        ::

             "sticky_left5"
             |
             |      "sticky_right5"
             |      |
            ---    ---
            GGGATCC
               TAGGTCA
               ----
                 |
                 "middle"



             "sticky_left3"
             |
             |      "sticky_right3"
             |      |
            ---    ---
               ATCCAGT
            CCCTAGG
               ----
                 |
                 "middle"



               "single_watson" (only an upper strand)
               |
            -------
            ATCCAGT
            |||||||



               "single_crick" (only a lower strand)
               |
            -------

            |||||||
            CCCTAGG


        Up to seven groups (0..6) are captured, but some are mutually exclusive
        which means that one of them is an empty string:

        0 or 1, not both, a DNA fragment has either 5' or 3' sticky end.

        2 or 5 or 6, a DNA molecule has a ds region or is single stranded.

        3 or 4, not both, either 5' or 3' sticky end.

        Note that internal single stranded regions are not identified and will
        be contained in the middle part if they are present.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> ds = Dseq("PPPATCFQZ")
        >>> ds
        Dseq(-9)
        GGGATC
           TAGTCA
        >>> parts = ds.get_parts()
        >>> parts
        DseqParts(sticky_left5='PPP', sticky_left3='', middle='ATC', sticky_right3='', sticky_right5='FQZ', single_watson='', single_crick='')
        >>> Dseq(parts.sticky_left5)
        Dseq(-3)
        GGG
        <BLANKLINE>
        >>> Dseq(parts.middle)
        Dseq(-3)
        ATC
        TAG
        >>> Dseq(parts.sticky_right5)
        Dseq(-3)
        <BLANKLINE>
        TCA

        Parameters
        ----------
        datastring : str
            A string with dscode.

        Returns
        -------
        namedtuple
            Seven string fields describing the DNA molecule.
            fragment(sticky_left5='', sticky_left3='',
                     middle='',
                     sticky_right3='', sticky_right5='',
                     single_watson='', single_crick='')

        """
        return get_parts(self._data.decode("ascii"))
