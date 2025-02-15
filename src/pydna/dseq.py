#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides the Dseq class for handling double stranded DNA sequences.

Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class
which can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology.
"""
import math as _math
import copy as _copy
import re as _re
import inspect as _inspect

from pydna.seq import Seq as _Seq

# from Bio.Seq import _translate_str, _SeqAbstractBaseClass
# from collections import namedtuple as _namedtuple
# from Bio.SeqFeature import SimpleLocation as _lc

from seguid import ldseguid as _ldseguid
from seguid import cdseguid as _cdseguid

from pydna._pretty import pretty_str as _pretty_str

# from pydna.utils import shift_location as _sl
from pydna.utils import rc as _rc
from pydna.utils import flatten as _flatten

from pydna.utils import cuts_overlap as _cuts_overlap
from pydna.utils import get_cutsite_pairs as _get_cutsite_pairs
from pydna.utils import complement as _complement
from pydna.utils import bp_dict
from pydna.utils import bp_dict_str
from pydna.utils import to_watson_table
from pydna.utils import to_crick_table
from pydna.utils import to_5tail_table
from pydna.utils import to_3tail_table
from pydna.utils import to_full_sequence
from pydna.utils import to_N as _to_N
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from pydna.common_sub_strings import terminal_overlap as _terminal_overlap

from Bio.Restriction import RestrictionBatch as _RestrictionBatch
from Bio.Restriction import CommOnly

from typing import (
    TYPE_CHECKING,
    List as _List,
    Tuple as _Tuple,
    Union as _Union,
    TypeVar as _TypeVar,
    Iterable as _Iterable,
)

try:
    from itertools import pairwise as _pairwise
except ImportError:

    def _pairwise(iterable):
        # pairwise('ABCDEFG') → AB BC CD DE EF FG
        iterator = iter(iterable)
        a = next(iterator, None)
        for b in iterator:
            yield a, b
            a = b


# pos_w_enz = _namedtuple("enzpos", "position enzyme recsite")
# nickpair = _namedtuple("nickpair", "position overhang enzyme recsite")

if TYPE_CHECKING:
    from Bio.Restriction import AbstractCut as _AbstractCut


# To represent any subclass of Dseq
DseqType = _TypeVar("DseqType", bound="Dseq")
EnzymesType = _TypeVar("EnzymesType", _RestrictionBatch, _Iterable["_AbstractCut"], "_AbstractCut")
CutSiteType = _Tuple[_Tuple[int, int], _Union["_AbstractCut", None]]


class Dseq(_Seq):
    """Dseq holds information for a double stranded (ds) DNA fragment that can be linear or circular.

    Dseq is a subclass of the Biopython Seq object. The DNA fragment can have single stranded (ss)
    regions, typically at either end. Such fragments are typically produced by restriction enzymes
    with staggered cuts.

    dsIUPAC [#]_ is an nn extension to the IUPAC alphabet used to describe ss regions:
    ::

            aaaGATC       GATCccc          ad-hoc representations
        CTAGttt               gggCTAG

        QFZJaaaPEXI       PEXIcccQFZJ      dsIUPAC



    Parameters
    ----------
    data : bytes, str
        a bytestring or string representing the ds DNA.

    circular : bool, optional
        True indicates that sequence is circular, False that it is linear.


    Examples
    --------
    The most common usage is probably to create a Dseq object as a
    part of a Dseqrecord object (see :class:`pydna.dseqrecord.Dseqrecord`).
    Dseqrecord object can hold metadata that closely follows the data stored in the
    Genbank sequence flat file format.

    There are three ways of creating a Dseq object directly listed below, but you can also
    use the function Dseq.from_full_sequence_and_overhangs() to create a Dseq:

    Only one argument (bytestring or string):

    >>> from pydna.dseq import Dseq
    >>> Dseq("aaa")
    Dseq(-3)
    aaa
    ttt

    The given string was interpreted as a blunt, linear double stranded
    nucleic acids fragment.

    We can use the dsIUPAC alphabet to create staggered DNA fragments:

    >>> from pydna.dseq import Dseq
    >>> Dseq("pexigggaaatqfzj")
    Dseq(-15)
    gatcgggaaat
        ccctttactag
    >>> Dseq("qfzjgggaaatpexi")
    Dseq(-15)
        gggaaatgatc
    ctagcccttta

    There is a check for internal consistency. A nuleic acid molecule has to have at
    least one phosphodiester bond in every position.
    ::

        GATT GATT
        CTAAG TAA

    >>> from pydna.dseq import Dseq
    >>> Dseq("GATTQPATT")
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/home/bjorn/python_packages/pydna/src/pydna/dseq.py", line 351, in __init__
        raise ValueError(f"Molecule is internally split: {', '.join([m.decode('ascii') for m in matches])}")
    ValueError: Molecule is internally split: QP

    A Dseq object can also be formed from the watson and crick strands using
    the pydna.dseq.pair helper function.

    The pair function accepts up to three arguments, watson, crick and ovhg=int.
    watson and crick are bytestrings or strings representing the watson (sense)
    or crick (antisense) strands of a double stranded DNA molecule.

    The ovhg parameter is a positive or negative integer describing the length
    of the crick strand overhang in at the left end (beginning) of the molecule.

    The ovhg parameter controls the stagger at the five prime end::

        dsDNA     ovhg (overhang)

          n...    2
        nnn...

         nn...    1
        nnn...

        nnn...    0
        nnn...

        nnn...   -1
         nn...

        nnn...   -2
          n...

    >>> from pydna.dseq import Dseq, pair
    >>> pair("GATT", "AATC", ovhg=0)
    b'GATT'
    >>> pair("aGATT", "AATC", ovhg=-1)
    b'eGATT'
    >>> pair("GATTa", "AATC", ovhg=0)
    b'GATTe'
    >>> pair("GATTa", "AATCc", ovhg=1)
    b'qGATTe'
    >>> Dseq(pair("GATTa", "AATCc", ovhg=1))
    Dseq(-6)
     GATTa
    cCTAA

    The ovhg parameter is optional. If both watson and crick are given, but not
    ovhg an attempt will be made to find the best annealing between the strands.
    There are limitations to this. For long fragments it is quite slow.

    The length of the annealing sequences have to be at least half the length of
    the shortest of the strands.

    >>> from pydna.dseq import Dseq, pair
    >>> pair("GATT", "AATC")
    b'GATT'
    >>> pair("aGATT", "AATC")
    b'eGATT'
    >>> pair("GATTa", "AATC")
    b'GATTe'
    >>> pair("GATTa", "AATCc")
    b'qGATTe'
    >>> Dseq(pair("GATTa", "AATCc"))
    Dseq(-6)
     GATTa
    cCTAA

    Example of creating Dseq objects with different amounts of stagger:

    >>> Dseq(pair(watson="agt", crick="actta", ovhg=-2))
    Dseq(-7)
    agt
      attca
    >>> Dseq(pair(watson="ata",crick="actta", ovhg=-1))
    Dseq(-6)
    ata
     attca
    >>> Dseq(pair(watson="taa",crick="actta",ovhg=0))
    Dseq(-5)
    taa
    attca
    >>> Dseq(pair(watson="aag",crick="actta",ovhg=1))
    Dseq(-5)
     aag
    attca
    >>> Dseq(pair(watson="agt",crick="actta",ovhg=2))
    Dseq(-5)
      agt
    attca

    The topology or shape of the fragment is set by circular = True or False

    >>> Dseq("aaa", circular = True)
    Dseq(o3)
    aaa
    ttt
    >>> Dseq("aaa", circular = False)
    Dseq(-3)
    aaa
    ttt

    >>> a=Dseq(pair("tttcccc","aaacccc"))
    >>> a
    Dseq(-11)
        tttcccc
    ccccaaa
    >>> a.ovhg
    4

    >>> b=Dseq(pair("ccccttt","ccccaaa"))
    >>> b
    Dseq(-11)
    ccccttt
        aaacccc
    >>> b.ovhg
    -4
    >>>

    Coercing to string

    >>> str(a.full_sequence)
    'ggggtttcccc'

    Coercing to string containing dsIUPAC code.

    >>> str(a)
    'qqqqtttiiii'

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
    >>> d=Dseq(s)
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


    The slicing of a circular Dseq object produce a linear Dseq object.


    >>> s="ggAtCc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggAtCc
    ccTaGg
    >>> d[1:5]
    Dseq(-4)
    gAtC
    cTaG


    The empty slice [:] returns a string or a linear Dseq unchanged, while
    a linear Dseq is returned for a circular Dseq object.
    This is the preferred way to linearize a circuler Dseq object.

    >>> s="ggatcc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggatcc
    cctagg
    >>> d[:]
    Dseq(-6)
    ggatcc
    cctagg
    >>>


    See Also
    --------
    pydna.dseqrecord.Dseqrecord
    .. [#] http://en.wikipedia.org/wiki/Klenow_fragment#The_exo-_Klenow_fragment

    """

    trunc = 30

    def __init__(
        self,
        data,
        circular=False,
        length=None,
        pos=0,
    ):
        super().__init__(data, length)
        t = self._data
        if circular:
            t += self._data[0:1]
        if matches := _re.findall(b"[PEXIpexi][QFZJqfzj]|[QFZJqfzj][PEXIpexi]", t):
            raise ValueError(f"Molecule is internally split: {', '.join([m.decode('ascii') for m in matches])}")
        self.circular = circular
        self.pos = pos

    @classmethod
    def quick(cls, data: bytes, *args, circular=False, pos=0, **kwargs):
        """Fastest way to instantiate an object of the Dseq class.

        No checks of parameters are made.
        """
        obj = cls.__new__(cls)  # Does not call Bio.Seq.__init__() which has lots of time consuming checks.
        obj.circular = circular
        obj.pos = pos
        obj._data = data
        return obj

    @classmethod
    def from_representation(cls, dsdna: str, *args, **kwargs):
        obj = cls.__new__(cls)
        obj.circular = False
        obj.pos = 0
        clean = _inspect.cleandoc("\n" + dsdna)
        watson, crick = [ln for ln in clean.splitlines() if ln.strip() and not ln.strip().startswith("Dseq(")]
        watson = f"{watson:<{len(crick)}}"
        crick = f"{crick:<{len(watson)}}"
        data = bytearray()
        for w, c in zip(watson, crick):
            try:
                data.extend(bp_dict[w.encode("ascii"), c.encode("ascii")])
            except KeyError as err:
                print(f"Base mismatch in representation {err}")
                raise
        obj._data = bytes(data)
        return obj

    @classmethod
    def from_full_sequence_and_overhangs(cls, full_sequence: str, crick_ovhg: int, watson_ovhg: int):
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
        watson = full_sequence
        crick = _complement(full_sequence)

        # If necessary, we trim the left side
        if crick_ovhg < 0:
            crick = -crick_ovhg * " " + crick[-crick_ovhg:]
        elif crick_ovhg > 0:
            watson = crick_ovhg * " " + watson[crick_ovhg:]

        # If necessary, we trim the right side
        if watson_ovhg < 0:
            watson = watson[:watson_ovhg] + " " * -watson_ovhg
        elif watson_ovhg > 0:
            crick = crick[:-watson_ovhg] + " " * watson_ovhg

        data = bytearray()

        for w, c in zip(watson, crick):
            data.extend(bp_dict[w.encode("ascii"), c.encode("ascii")])

        return cls(data)

    def is_watson(self):
        return not any(c in "QFZJqfzjGATCOUgatcou" for c in self._data.decode("ascii"))

    def is_crick(self):
        return not any(c in "PEXIpexiGATCOUgatcou" for c in self._data.decode("ascii"))

    def is_ds(self):
        return not (self.is_watson() or self.is_crick())

    def mw(self) -> float:
        """This method returns the molecular weight of the DNA molecule
        in g/mol. The following formula is used::

               MW = (A x 313.2) + (T x 304.2) +
                    (C x 289.2) + (G x 329.2) +
                    (N x 308.9) + 79.0
        """
        nts = f"{self.watson}{self.crick}".lower()
        return (
            313.2 * nts.count("a")
            + 304.2 * nts.count("t")
            + 289.2 * nts.count("c")
            + 329.2 * nts.count("g")
            + 308.9 * nts.count("n")
            + 79.0
        )

    @property
    def watson(self):
        return self._data.translate(to_watson_table).strip().decode("ascii")

    @property
    def crick(self):
        return self._data.translate(to_crick_table).strip().decode("ascii")[::-1]

    @property
    def ovhg(self):
        ohw = len(_re.match(b"^[PEXIpexi]*", self._data).group(0))
        ohc = len(_re.match(b"^[QFZJqfzj]*", self._data).group(0))
        return -ohw or ohc

    @property  # TODO: talk about the naming
    def full_sequence(self):
        return self.__class__(self._data.translate(to_full_sequence), circular=self.circular)

    def watson_ovhg(self) -> int:
        """Returns the overhang of the watson strand at the three prime."""
        ohw = len(_re.search(b"[PEXIpexi]*$", self._data).group(0))
        ohc = len(_re.search(b"[QFZJqfzj]*$", self._data).group(0))
        return ohw or -ohc

    def right_end_position(self) -> _Tuple[int, int]:
        """The index in the full sequence of the watson and crick end positions.

        full sequence (str(self)) for all three cases is AAA

        ```
        AAA               AA                   AAA
        TT                TTT                  TTT
        Returns (3, 2)    Returns (2, 3)       Returns (3, 3)
        ```

        """
        if self.watson_ovhg() < 0:
            return len(self) + self.watson_ovhg(), len(self)
        return len(self), len(self) - self.watson_ovhg()

    def left_end_position(self) -> _Tuple[int, int]:
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

    def five_prime_end(self) -> _Tuple[str, str]:
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.five_prime_end()
        ('blunt', b'')
        >>> b=Dseq("QZJaaa")
        >>> b
        Dseq(-6)
           aaa
        CAGttt
        >>> b.five_prime_end()
        ("3'", b'gac')
        >>> c=Dseq("PXEaaa")
        >>> c
        Dseq(-6)
        GTAaaa
           ttt
        >>> c.five_prime_end()
        ("5'", b'gta')
        >>>

        See also
        --------
        pydna.dseq.Dseq.three_prime_end

        """
        sticky5 = _re.match(b"^[PEXIpexi]*", self._data).group(0).translate(to_watson_table).lower()
        sticky3 = _re.match(b"^[QFZJqfzj]*", self._data).group(0).translate(to_crick_table).lower()

        type_ = "blunt"

        assert (sticky5 or sticky3) or (not sticky5 and not sticky3), "End has to be either 5' or 3', not both."

        if sticky5:
            type_ = "5'"

        if sticky3:
            type_ = "3'"

        return type_, sticky5 or sticky3[::-1]

    def three_prime_end(self) -> _Tuple[str, str]:
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', b'')
        >>> b=Dseq("aaaQZJ")
        >>> b
        Dseq(-6)
        aaa
        tttCAG
        >>> b.three_prime_end()
        ("5'", b'gac')
        >>> c=Dseq("aaaPXE")
        >>> c
        Dseq(-6)
        aaaGTA
        ttt
        >>> c.three_prime_end()
        ("3'", b'gta')
        >>>

        See also
        --------
        pydna.dseq.Dseq.five_prime_end

        """

        sticky3 = _re.search(b"[PEXIpexi]*$", self._data).group(0).translate(to_watson_table).lower()
        sticky5 = _re.search(b"[QFZJqfzj]*$", self._data).group(0).translate(to_crick_table).lower()

        type_ = "blunt"

        assert (sticky5 or sticky3) or (not sticky5 and not sticky3), "End has to be either 5' or 3', not both."

        if sticky5:
            type_ = "5'"

        if sticky3:
            type_ = "3'"

        return type_, sticky5[::-1] or sticky3

    def translate(self, *args, **kwargs):  # TODO: discuss this
        obj = _Seq(self.full_sequence._data)
        return obj.translate(*args, **kwargs)

    def __str__(self):
        return self._data.decode("ascii")  # self.full_sequence._data.decode("ascii")

    # def __getitem__(self, sl: slice) -> "Dseq":
    #     """Returns a subsequence. This method is used by the slice notation"""
    #     if isinstance(sl, slice):
    #         if sl.start == sl.stop and self.circular and (sl.start or sl.start==0) and (sl.stop or sl.stop==0):
    #             return self[sl.start:] + self[:sl.start]
    #     return self.__class__(super().__getitem__(sl))

    def __getitem__(self, sl: slice) -> "Dseq":
        """Returns a subsequence. This method is used by the slice notation"""
        if isinstance(sl, slice):
            if sl.start is None and sl.stop is None and self.circular:
                return self.__class__(self._data[sl])
        return self.__class__(super().__getitem__(sl))

    def __contains__(self, item):
        item = item if isinstance(item, bytes) else item.encode("ascii")
        return self.full_sequence._data.__contains__(item)

    def __hash__(self):
        return hash((self._data, self.circular))

    def __eq__(self, other: DseqType) -> bool:
        """Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties. This comparison is case
        insensitive.

        """
        try:
            same = self.__hash__() == other.__hash__()
            # same = (
            #     other.watson.lower() == self.watson.lower()
            #     and other.crick.lower() == self.crick.lower()
            #     and other.ovhg == self.ovhg
            #     and self.circular == other.circular
            # )
            # Also test for alphabet ?
        except AttributeError:
            same = False
        return same

    def __repr__(self):
        header = f"{self.__class__.__name__}({({False: '-', True: 'o'}[self.circular])}{len(self)})"
        m = _re.match(
            b"([PEXIpexi]*)([QFZJqfzj]*)(?=[GATCUOgatcuo])(.*)(?<=[GATCUOgatcuo])([PEXIpexi]*)([QFZJqfzj]*)|([PEXIpexiQFZJqfzj]+)",
            self._data,
        )
        result = m.groups() if m else (b"",) * 6
        sticky_left5, sticky_left3, middle, sticky_right5, sticky_right3, single = result
        if len(self) > self.trunc:
            sticky_left5 = (
                sticky_left5[:4] + b"22" + sticky_left5[-4:]
                if sticky_left5 and len(sticky_left5) > 10
                else sticky_left5
            )
            sticky_left3 = (
                sticky_left3[:4] + b"11" + sticky_left3[-4:]
                if sticky_left3 and len(sticky_left3) > 10
                else sticky_left3
            )
            middle = middle[:4] + b".." + middle[-4:] if middle and len(middle) > 10 else middle
            sticky_right5 = (
                sticky_right5[:4] + b"22" + sticky_right5[-4:]
                if sticky_right5 and len(sticky_right5) > 10
                else sticky_right5
            )
            sticky_right3 = (
                sticky_right3[:4] + b"11" + sticky_right3[-4:]
                if sticky_right3 and len(sticky_right3) > 10
                else sticky_right3
            )
        r = (sticky_left5 or sticky_left3 or b"") + (middle or b"") + (sticky_right5 or sticky_right3 or single or b"")
        return _pretty_str(
            f"{header}\n{r.translate(to_watson_table).decode().rstrip()}\n{r.translate(to_crick_table).decode().rstrip()}"
        )

        # return _pretty_str(
        #     f"{header}\n{r.translate(to_watson_table).decode().rstrip()}\n{_complement(r.translate(to_crick_table)).decode().rstrip()}"
        # )

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
        return Dseq(_rc(self._data), circular=self.circular)

    rc = reverse_complement  # alias for reverse_complement

    def shifted(self: DseqType, shift: int) -> DseqType:
        """Shifted version of a circular Dseq object."""
        if not self.circular:
            raise TypeError("DNA is not circular.")
        shift = shift % len(self)
        if not shift:
            return _copy.deepcopy(self)
        else:
            return self.__class__(self._data[shift:] + self._data[:shift], circular=True)

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
        >>> c = Dseq("ietcgazj")
        >>> c
        Dseq(-8)
        catcga
          agctag
        >>> c.looped()
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/usr/local/lib/python2.7/dist-packages/pydna/dsdna.py", line 357, in looped
            if type5 == type3 and str(sticky5) == str(rc(sticky3)):
        TypeError: DNA cannot be circularized.
        5' and 3' sticky ends not compatible!
        >>>

        """
        if self.circular:
            return _copy.deepcopy(self)

        ms = _re.fullmatch(b"(?:.*)([GATCgatc])([QFZJqfzj]+|[PEXIpexi]+)", self._data)
        mo = _re.fullmatch(b"([PEXIpexi]+|[QFZJqfzj]+)([GATCgatc])(?:.*)", self._data)

        sticky_left = ms.group(2) if ms else b""
        sticky_right = mo.group(1) if mo else b""

        sticky_left_just = bytes(f"{sticky_left.decode():<{len(sticky_right)}}", "utf8")
        sticky_right_just = bytes(f"{sticky_right.decode():>{len(sticky_left)}}", "utf8")

        assert len(sticky_left_just) == len(sticky_right_just)

        junction = b"".join(
            [bp_dict.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_left_just, sticky_right_just)]
        )

        if b"-" in junction:
            raise TypeError("DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!")

        return self.__class__(
            junction + self._data[len(sticky_left) or None : -len(sticky_right) or None], circular=True
        )

    def tolinear(self: DseqType) -> DseqType:  # pragma: no cover
        """Returns a blunt, linear copy of a circular Dseq object. This can
        only be done if the Dseq object is circular, otherwise a
        TypeError is raised.

        This method is deprecated, use slicing instead. See example below.

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc", circular=True)
        >>> a
        Dseq(o8)
        catcgatc
        gtagctag
        >>> a[:]
        Dseq(-8)
        catcgatc
        gtagctag
        >>>

        """
        import warnings as _warnings
        from pydna import _PydnaDeprecationWarning

        _warnings.warn(
            "tolinear method is obsolete; " "please use obj[:] " "instead of obj.tolinear().",
            _PydnaDeprecationWarning,
        )
        if not self.circular:
            raise TypeError("DNA is not circular.\n")
        selfcopy = _copy.deepcopy(self)
        selfcopy.circular = False
        return selfcopy

    # def _add(self: DseqType, other: DseqType, perfectmatch) -> DseqType:

    #     if self.circular:
    #         raise TypeError("circular DNA cannot be ligated!")
    #     try:
    #         if other.circular:
    #             raise TypeError("circular DNA cannot be ligated!")
    #     except AttributeError:
    #         pass

    #     if not other:
    #         return _copy.deepcopy(self)
    #     elif not self:
    #         return _copy.deepcopy(other)

    #     ms = _re.fullmatch(b"(?:.*)([GATCgatc])([QFZJqfzj]+|[PEXIpexi]+)", self._data)
    #     mo = _re.fullmatch(b"([PEXIpexi]+|[QFZJqfzj]+)([GATCgatc])(?:.*)", other._data)

    #     # ms = _re.search(b"([PEXIpexiQFZJqfzj])$", self._data)
    #     # mo = _re.match(b"^([PEXIpexiQFZJqfzj]+)", other._data)

    #     sticky_self = ms.group(2) if ms else b""
    #     sticky_other = mo.group(1) if mo else b""

    #     sticky_self_just = bytes(f"{sticky_self.decode():<{len(sticky_other)}}", 'utf8')
    #     sticky_other_just = bytes(f"{sticky_other.decode():>{len(sticky_self)}}", 'utf8')

    #     assert len(sticky_self_just) == len(sticky_other_just)

    #     junction = b"".join(
    #         [bp_dict.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_self_just, sticky_other_just)]
    #     )

    #     if b"-" in junction:
    #         raise TypeError("sticky ends not compatible!")

    #     return self.__class__(
    #         self._data[: -len(sticky_self) or None] + junction + other._data[len(sticky_other) or None :]
    #     )

    def _add(self: DseqType, other: DseqType, perfectmatch) -> DseqType:

        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        if not other:
            return _copy.deepcopy(self)
        elif not self:
            return _copy.deepcopy(other)

        ms = _re.search(b"([PEXIpexiQFZJqfzj]+)$", self._data)
        mo = _re.match(b"^([PEXIpexiQFZJqfzj]+)", other._data)

        sticky_self = ms.group() if ms else b""
        sticky_other = mo.group() if mo else b""

        sticky_self_just = bytes(f"{sticky_self.decode():<{len(sticky_other)}}", "utf8")
        sticky_other_just = bytes(f"{sticky_other.decode():>{len(sticky_self)}}", "utf8")

        assert len(sticky_self_just) == len(sticky_other_just)

        if perfectmatch:

            junction = b"".join(
                [bp_dict.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_self_just, sticky_other_just)]
            )

            if b"-" in junction:
                raise TypeError("sticky ends not compatible!")

        else:

            result = _terminal_overlap(
                sticky_self.translate(to_full_sequence).decode("ascii"),
                sticky_other.translate(to_full_sequence).decode("ascii") + "@",
                limit=1,
            )

            if not result:
                return None

            (startx, starty, length), *_ = result

            sticky_self = sticky_self[startx : startx + length]
            sticky_other = sticky_other[starty : starty + length]

            junction = b"".join(
                [bp_dict.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_self, sticky_other)]
            )

        return self.__class__(
            self._data[: -len(sticky_self) or None] + junction + other._data[len(sticky_other) or None :]
        )

    def __add__(self: DseqType, other: DseqType) -> DseqType:
        return self._add(other, perfectmatch=True)

    def __truediv__(self, other):
        # x / y   	__rfloordiv__  https://www.pythonmorsels.com/every-dunder-method
        return self._add(other, perfectmatch=False)

    def __invert__(self):
        # ~x https://www.pythonmorsels.com/every-dunder-method
        return self.reverse_complement()

    def __mul__(self: DseqType, number: int) -> DseqType:
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(number)))
        if number <= 0:
            return self.__class__("")
        new = _copy.deepcopy(self)
        for i in range(number - 1):
            new += self
        return new

    def _fill_in_three_prime(self: DseqType, nucleotides: str) -> str:
        """
        PEXINNNN

        GATCNNNN
            NNNN

        GATCNNNN
        CTAGNNNN

        """

        nucleotides = nucleotides if isinstance(nucleotides, bytes) else nucleotides.encode("ascii")

        if not nucleotides:
            return _copy.deepcopy(self)

        s = _complement(nucleotides.upper() + nucleotides.lower()).translate(to_3tail_table)

        def repl(m):
            return m.group(1).translate(to_watson_table)

        return self.__class__(_re.sub(b"([%b]+)(?=[GATCgatc])" % s, repl, self._data), circular=False)

    def _fill_in_five_prime(self: DseqType, nucleotides: str) -> str:
        """
        NNNNQFZJ

        NNNN----
        NNNNCTAG

        NNNNGATC
        NNNNCTAG

        """
        nucleotides = nucleotides if isinstance(nucleotides, bytes) else nucleotides.encode("ascii")

        if not nucleotides:
            return _copy.deepcopy(self)

        s = (nucleotides.upper() + nucleotides.lower()).translate(to_5tail_table)

        def repl(m):
            return m.group(1).translate(to_full_sequence)

        return self.__class__(_re.sub(b"(?<=[GATCgatc])([%b]+)" % s, repl, self._data), circular=False)

    def fill_in(self, nucleotides: _Union[bytes, str]) -> "Dseq":
        """Fill in of five prime protruding end with a DNA polymerase
        that has only DNA polymerase activity (such as exo-klenow [#]_)
        and any combination of A, G, C or T.

        Parameters
        ----------

        nucleotides : str

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.fill_in("GATC")
        Dseq(-3)
        aaa
        ttt
        >>> b = Dseq("iaaaq")
        >>> b
        Dseq(-5)
        caaa
         tttc
        >>> b.fill_in("g")
        Dseq(-5)
        caaag
        gtttc
        >>> b.fill_in("tac")
        Dseq(-5)
        caaa
         tttc
        >>> c=Dseq("jaaai")
        >>> c
        Dseq(-5)
         aaac
        gttt
        >>> c.fill_in("gatc")  # DNA pol does not process 3'-5'
        Dseq(-5)
         aaac
        gttt
        >>>

        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Klenow_fragment#The_exo-_Klenow_fragment

        """

        return self._fill_in_five_prime(nucleotides)._fill_in_three_prime(nucleotides)

    def mung(self) -> "Dseq":
        """
        Simulates treatment a nuclease with 5'-3' and 3'-5' single
        strand specific exonuclease activity (such as mung bean nuclease [#]_)

        ::

             ggatcc    ->     gatcc
              ctaggg          ctagg

              ggatcc   ->      ggatc
             tcctag            cctag

         >>> from pydna.dseq import Dseq
         >>> b=Dseq("iaaaq")
         >>> b
         Dseq(-5)
         caaa
          tttc
         >>> b.mung()
         Dseq(-3)
         aaa
         ttt
         >>> c=Dseq("jaaai")
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

        return self.__class__(self._data.strip(b"PEXIQFZJpexiqfzj"), circular=False)

    def T4(self, nucleotides: _Union[bytes, str]) -> "Dseq":
        """Fill in five prime protruding ends and chewing back
        three prime protruding ends by a DNA polymerase providing both
        5'-3' DNA polymerase activity and 3'-5' nuclease acitivty
        (such as T4 DNA polymerase). This can be done in presence of any
        combination of the four A, G, C or T. Removing one or more nucleotides
        can facilitate engineering of sticky ends. Default are all four nucleotides together.

        Parameters
        ----------
        nucleotides : str


        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("gatcgatc")
        >>> a
        Dseq(-8)
        gatcgatc
        ctagctag
        >>> a.T4("t")
        Dseq(-8)
        gatcgat
         tagctag
        >>> a.T4("a")
        Dseq(-8)
        gatcga
          agctag
        >>> a.T4("g")
        Dseq(-8)
        gatcg
           gctag
        >>>

        """

        nucleotides = nucleotides if isinstance(nucleotides, bytes) else nucleotides.encode("ascii")
        to_remove = bytes(set(b"GATCgatc") - set(nucleotides.upper() + nucleotides.lower()))
        strp = self.lstrip(b"QFZJqfzj").rstrip(b"PEXIpexi")  # remove 3' sticky ends on both sides

        if not to_remove:
            return strp.fill_in(nucleotides)

        def repl1(m):
            return m.group(1).translate(to_3tail_table)

        new = _re.sub(b"(?:(?<=[PEXIpexi])|^)([%b]+)" % _complement(to_remove), repl1, strp._data)

        def repl2(m):
            return m.group(1).translate(to_5tail_table)

        final = _re.sub(b"([%b]+)(?=([QFZJqfzj]|$))" % to_remove, repl2, new)

        if not _re.search(b"([PEXIpexi]|[QFZJqfzj])[GATCgatc]+([PEXIpexi]|[QFZJqfzj])", final):
            import warnings as _warnings
            from pydna import _PydnaWarning

            _warnings.warn(
                "DNA fragment was completely digested by the 3'-5' nuclease activity of T4 DNA pol", _PydnaWarning
            )
            final = b""

        return self.__class__(final).fill_in(nucleotides)

    t4 = T4  # alias for the T4 method.

    def exo1_front(self: DseqType, n=1) -> DseqType:
        """5'-3' resection at the start (left side) of the molecule."""
        d = _copy.deepcopy(self)
        d.ovhg += n
        d.watson = d.watson[n:]
        return d

    def exo1_end(self: DseqType, n=1) -> DseqType:
        """5'-3' resection at the end (right side) of the molecule."""
        d = _copy.deepcopy(self)
        d.crick = d.crick[n:]
        return d

    def no_cutters(self, batch: _Union[_RestrictionBatch, None] = None) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch not cutting sequence."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if not sitelist}
        return _RestrictionBatch(ncut)

    def unique_cutters(self, batch: _Union[_RestrictionBatch, None] = None) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence once."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=1, batch=batch)

    once_cutters = unique_cutters  # alias for unique_cutters

    def twice_cutters(self, batch: _Union[_RestrictionBatch, None] = None) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence twice."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=2, batch=batch)

    def n_cutters(self, n=3, batch: _Union[_RestrictionBatch, None] = None) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting n times."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if len(sitelist) == n}
        return _RestrictionBatch(ncut)

    def cutters(self, batch: _Union[_RestrictionBatch, None] = None) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence at least once."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if sitelist}
        return _RestrictionBatch(ncut)

    def seguid(self) -> str:
        """SEGUID checksum for the sequence."""
        if self.circular:
            cs = _cdseguid(self.watson.upper(), self.crick.upper(), alphabet="{DNA-extended},AU")
        else:
            """docstring."""
            oh_L = self.ovhg
            oh_R = self.watson_ovhg()
            w = oh_L * "-" + self.watson.upper().replace(" ", "-") + -oh_R * "-"
            c = oh_R * "-" + self.crick.upper().replace(" ", "-") + -oh_L * "-"
            cs = _ldseguid(w, c, alphabet="{DNA-extended},AU")
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
        >>> a=Dseq("jgat")
        >>> a
        Dseq(-4)
         gat
        gcta
        >>> a.isblunt()
        False
        >>> a=Dseq("gatj")
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
        return (
            (self._data[0] not in b"PEXIQFZJpexiqfzj")
            and (self._data[-1] not in b"PEXIQFZJpexiqfzj")
            and not self.circular
        )

    def cas9(self, RNA: str) -> _Tuple[slice, ...]:
        """docstring."""
        bRNA = bytes(RNA, "ASCII")
        slices = []
        cuts = [0]
        for m in _re.finditer(bRNA, self._data):
            cuts.append(m.start() + 17)
        cuts.append(len(self))
        slices = tuple(slice(x, y, 1) for x, y in zip(cuts, cuts[1:]))
        return slices

    def terminal_transferase(self, nucleotide="a") -> "Dseq":
        """docstring."""
        first = _complement(nucleotide).translate(to_5tail_table)
        last = nucleotide.translate(to_3tail_table)
        return self.__class__(first + self._data + last)

    def user(self):
        """
        USER Enzyme is a mixture of Uracil DNA glycosylase (UDG) and the
        DNA glycosylase-lyase Endonuclease VIII. UDG catalyses the excision
        of an uracil base, forming an abasic (apyrimidinic) site while leaving
        the phosphodiester backbone intact (2,3).

        Returns
        -------
        None.

        """
        return Dseq(self._data.translate(bytes.maketrans(b"UuOo", b"FfEe")))

    def cut(self: DseqType, *enzymes: EnzymesType) -> _Tuple[DseqType, ...]:
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

        if not cutsites:
            return tuple()

        cutsite_pairs, shift, stuffer = _get_cutsite_pairs(cutsites, self.circular, len(self))

        return tuple(self.apply_cut(*cutsite_pair, shift, stuffer) for cutsite_pair in cutsite_pairs)

    def melt(self, length):
        if not length or length < 1:
            return tuple()

        new, strands = self.shed_ss_dna(length)

        cutsites = new.get_ds_meltsites(length)

        cutsite_pairs, shift, stuffer = _get_cutsite_pairs(cutsites, self.circular, len(self))

        result = tuple(new.apply_cut(*cutsite_pair, shift, stuffer) for cutsite_pair in cutsite_pairs)

        result = tuple([new]) if strands and not result else result

        return strands + result

    def shed_ss_dna(self, length):

        new, strands, intervals = self._shed_ss_dna(length)

        return Dseq(new), strands

    def _shed_ss_dna(self, length):

        watsonnicks, cricknicks = self.get_ss_meltsites(length)

        watsonstrands = []
        crickstrands = []

        new = self._data

        for x, y in watsonnicks:
            stuffer = new[x:y]
            ss = Dseq.quick(stuffer.translate(to_5tail_table))
            new = new[:x] + stuffer.translate(to_3tail_table) + new[y:]
            watsonstrands.append((x, y, ss))

        for x, y in cricknicks:
            stuffer = new[x:y]
            ss = Dseq.quick(stuffer.translate(to_3tail_table))
            new = new[:x] + stuffer.translate(to_5tail_table) + new[y:]
            crickstrands.append((x, y, ss))

        ordered_strands = sorted(watsonstrands + crickstrands)

        strands = []

        for x, y, ss in ordered_strands:
            seq = ss._data[::-1].translate(to_watson_table).strip() or ss._data.translate(to_crick_table).strip()
            strands.append(_Seq(seq))

        intervals = tuple((x, y) for x, y, ss in ordered_strands)

        return Dseq(new), tuple(strands), intervals

    def apply_cut(self, left_cut: CutSiteType, right_cut: CutSiteType, shift, stuffer) -> "Dseq":
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
        >>> from pydna.utils import get_cutsite_pairs
        >>> dseq = Dseq('aaGAATTCaaGAATTCaa')
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((3, -4), EcoRI), ((11, -4), EcoRI)]
        >>> pairs, shift, stuffer = get_cutsite_pairs(cutsites, circular=False, length=len(dseq))
        >>> p1, p2, p3 = pairs
        >>> p1
        (((0, 0), None), ((3, -4), EcoRI))
        >>> dseq.apply_cut(*p1, shift, stuffer)
        Dseq(-7)
        aaG
        ttCTTAA
        >>> p2
        (((3, -4), EcoRI), ((11, -4), EcoRI))
        >>> dseq.apply_cut(*p2, shift, stuffer)
        Dseq(-12)
        AATTCaaG
            GttCTTAA
        >>> p3
        (((11, -4), EcoRI), ((18, 0), None))
        >>> dseq.apply_cut(*p3, shift, stuffer)
        Dseq(-7)
        AATTCaa
            Gtt

        >>> dseq = Dseq('TTCaaGAA', circular=True)
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((6, -4), EcoRI)]
        >>> pair, shift, stuffer = get_cutsite_pairs(cutsites, circular=True, length=len(dseq))
        >>> pair
        [(((0, -4), EcoRI), ((8, -4), EcoRI))]
        >>> dseq.apply_cut(*pair[0], shift, stuffer)
        Dseq(-12)
        AATTCaaG
            GttCTTAA

        """
        if _cuts_overlap(left_cut, right_cut, len(self), self.circular):
            raise ValueError("Cuts by {} {} overlap.".format(left_cut[1], right_cut[1]))

        (left_watson_cut, left_overhang), left_meta_data = left_cut
        (right_watson_cut, right_overhang), right_meta_data = right_cut

        table1 = to_5tail_table if left_overhang > 0 else to_3tail_table
        table2 = to_5tail_table if right_overhang < 0 else to_3tail_table

        left_stck_begin = min(left_watson_cut, left_watson_cut - left_overhang)
        left_stck_end = left_stck_begin + abs(left_overhang)

        right_stck_begin = min(right_watson_cut, right_watson_cut - right_overhang)
        right_stck_end = right_stck_begin + abs(right_overhang)

        cutfrom = self.shifted(shift)[:] if shift else self[:]
        cutfrom = (cutfrom + cutfrom[:stuffer])._data

        left_sticky_end = cutfrom[left_stck_begin:left_stck_end].translate(table1)
        ds_middle_part = cutfrom[left_stck_end:right_stck_begin]
        right_sticky_end = cutfrom[right_stck_begin:right_stck_end].translate(table2)

        return Dseq(left_sticky_end + ds_middle_part + right_sticky_end)

    def get_cutsites(self: DseqType, *enzymes: EnzymesType) -> _List[CutSiteType]:
        """Enzyme cutsites

        Returns a sorted list (by distance from the left) of tuples. Each tuple (`((cut_watson, ovhg), enz)`)
        contains a tuple containing the cut position and the overhang value for the enzyme followed by the
        an enzyme object.

        - `cut_watson` is a positive integer contained in `[0, len(seq))`, where `seq` is the sequence
          that will be cut. It represents the position of the cut on the watson strand, using the full
          sequence as a reference. By "full sequence" I mean the one you would get from `str(Dseq)`.
        - `ovhg` is the overhang left after the cut. It has the same meaning as `ovhg` in
          the `Bio.Restriction` enzyme objects, or pydna's `Dseq` property (see comments in code below).
        - `enz` is the enzyme object. It's not necessary to perform the cut, but can be
           used to keep track of which enzyme was used.

        Cuts are only returned if
        - Both cut positions fall inside the sequence (could be moved to Biopython).
        - Both cut positions fall double stranded region od the sequence.
        - At least one cut pair fulfils the criteria above for enzymes that cut twice (such as BcgI).

        Parameters
        ----------

        enzymes : Union[_RestrictionBatch,list[_AbstractCut]]

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

        if len(enzymes) == 1 and isinstance(enzymes[0], _RestrictionBatch):
            # argument is probably a RestrictionBatch
            enzymes = [e for e in enzymes[0]]

        enzymes = _flatten(enzymes)

        # The out list will eventually contain all valid cuts.
        out = []

        # The restriction enzyme search method only allows letters in the IUPAC extended alphabet.
        # See the Bio.Restriction.FormattedSequence for details.
        # For this reason we need to construct a temporary sequence where we filter and rplace with N
        model = Dseq.quick(self._data.translate(_to_N))
        # The ln variable contains the length so that we do not have to calculate it repeatedly.
        ln = len(self)

        for e in enzymes:

            # The ovhg property of the enzyme indicate the offset between the
            # cut in the watson strand and the cut in the complementary crick
            # strand. This can be positive or negative as indicated below for
            # three common enzymes that cut the same sequence in different ways.

            oh = e.ovhg

            # G|G T A C C    Acc65I.ovhg = -4
            #  ---------
            # C C A T G|G

            # G G T|A C C    NlaIV.ovhg = 0
            # C C A|T G G

            # G G T A C|C    KpnI.ovhg = 4
            #  ---------
            # C|C A T G G

            # Positions of the cut on the watson strand. They are 1-based, so we subtract
            # 1 to get 0-based positions
            cuts_watson = [c - 1 for c in e.search(model, linear=not self.circular)]
            # The search method of the enzyme can handle linear or circular sequences
            # Cuts in circular sequences are reported even enzyme straddles the sequence origin

            if self.circular:
                # If the sequence is circular, we need to construct a temporary sequence
                # with a stuffer on each end in case the enzyme overhang cuts before or after
                # the origin
                spacer = abs(oh)
                cutfrom = Dseq.quick(self._data[-spacer:] + self._data + self._data[:spacer])
                cuts_watson_inside = cuts_watson
            else:
                # If the sequence is linear, we need to filter cuts that cut outside
                spacer = 0
                cuts_watson_inside = [c for c in cuts_watson if 0 <= c - oh <= ln]
                # Here we cut from original sequence, no spacers needed
                cutfrom = self

            for cut_on_watson in cuts_watson_inside:
                begin, end = sorted((cut_on_watson, cut_on_watson - oh))
                # We extract the sticky end from the temporary sequence
                sticky_end = cutfrom[begin + spacer : end + spacer]
                # If sticky end is all double stranded, the cut is valid.
                # This implies that the sequence does not contain any
                # letter associated with single stranded DNA (dsIUPAC)
                if all(c not in "PEXIpexiQFZJqfzj" for c in sticky_end):
                    out.append(((cut_on_watson, oh), e))

        return sorted(out)

    def get_ss_meltsites(self: DseqType, length) -> _List[CutSiteType]:

        if length < 1:
            return ()

        regex = _re.compile(
            (
                f"(?P<watson>((?<=[PEXIpexi]))([GATCgatcUuOo]{{1,{length}}})((?=[^QFZJqfzjGATCgatcUuOo])))|"
                f"(?P<crick>((?<=[QFZJqfzj]))([GATCgatcUuOo]{{1,{length}}})((?=[^PEXIpexiGATCgatcUuOo])))"
            ).encode("ascii")
        )

        if self.circular:

            spacer = length

            cutfrom = self._data[-length:] + self._data + self._data[:length]

        else:

            spacer = 0

            cutfrom = self._data

        watsoncuts = []
        crickcuts = []

        for m in regex.finditer(cutfrom):

            if m.lastgroup == "watson":
                cut1 = m.start() + spacer
                cut2 = m.end() + spacer
                watsoncuts.append((cut1, cut2))
            else:
                assert m.lastgroup == "crick"
                cut1 = m.start() + spacer
                cut2 = m.end() + spacer
                crickcuts.append((cut1, cut2))

        return watsoncuts, crickcuts

    def get_ds_meltsites(self: DseqType, length) -> _List[CutSiteType]:
        """Double stranded DNA melt sites

        Returns a sorted list (by distance from the left) of tuples. Each tuple (`((cut_watson, ovhg), enz)`)
        contains a tuple containing the cut position and the overhang value for the enzyme followed by the
        an enzyme object.

        """

        if length < 1:
            return ()

        regex = _re.compile(
            (
                f"(?P<watson>((?<=[PEXIpexi])|^)([GATCgatcUuOo]{{1,{length}}})((?=[^PEXIpexiGATCgatcUuOo])|$))|"
                f"(?P<crick>((?<=[QFZJqfzj])|^)([GATCgatcUuOo]{{1,{length}}})((?=[^QFZJqfzjGATCgatcUuOo])|$))"
            ).encode("ascii")
        )

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

    # def cutsite_is_valid(self, cutsite: CutSiteType) -> bool:
    #     """Returns False if:
    #     - Cut positions fall outside the sequence (could be moved to Biopython)
    #     - Overhang is not double stranded
    #     - Recognition site is not double stranded or is outside the sequence
    #     - For enzymes that cut twice, it checks that at least one possibility is valid
    #     """

    #     assert cutsite is not None, "cutsite is None"

    #     enz = cutsite[1]
    #     watson, crick, ovhg = self.get_cut_parameters(cutsite, True)

    #     # The overhang is double stranded
    #     overhang_dseq = self[watson:crick] if ovhg < 0 else self[crick:watson]
    #     if overhang_dseq.ovhg != 0 or overhang_dseq.watson_ovhg() != 0:
    #         return False

    #     # The recognition site is double stranded and within the sequence
    #     start_of_recognition_site = watson - enz.fst5
    #     if start_of_recognition_site < 0:
    #         start_of_recognition_site += len(self)
    #     end_of_recognition_site = start_of_recognition_site + enz.size
    #     if self.circular:
    #         end_of_recognition_site %= len(self)
    #     recognition_site = self[start_of_recognition_site:end_of_recognition_site]
    #     if len(recognition_site) == 0 or recognition_site.ovhg != 0 or recognition_site.watson_ovhg() != 0:
    #         if enz is None or enz.scd5 is None:
    #             return False
    #         else:
    #             # For enzymes that cut twice, this might be referring to the second one
    #             start_of_recognition_site = watson - enz.scd5
    #             if start_of_recognition_site < 0:
    #                 start_of_recognition_site += len(self)
    #             end_of_recognition_site = start_of_recognition_site + enz.size
    #             if self.circular:
    #                 end_of_recognition_site %= len(self)
    #             recognition_site = self[start_of_recognition_site:end_of_recognition_site]

    #             if len(recognition_site) == 0 or recognition_site.ovhg != 0 or recognition_site.watson_ovhg() != 0:
    #                 return False

    #     return True


def pair(watson: [str, bytes], crick: [str, bytes], ovhg: int = None):
    """
    Pair two DNA strands, Watson and Crick, by finding their longest overlapping sequence.

    The function returns the paired (annealed) DNA sequence in dsIUPAC as an ASCII encoded
    bytes object.

    Parameters
    ----------
    watson : str or bytes
        The first DNA strand (Watson). Can be a string or a bytes object.
        If bytes, it will be decoded to ASCII.
    crick : str or bytes
        The second DNA strand (Crick). Can be a string or a bytes object.
        If bytes, it will be decoded to ASCII.
    ovhg : int
        Amount of stagger

    Returns
    -------
    bytes
        A bytes object representing the paired DNA sequence, encoded in ASCII.

    Raises
    ------
    ValueError
        If no overlapping sequences are found, or if there are multiple possible alignments
        of the same maximum overlap length.

    Notes
    -----
    - The function determines the overlap between the two DNA strands by identifying
      the longest common substring between the lowercase Watson strand and the reverse
      complement of the Crick strand.
    - If multiple overlaps of the same maximum length exist, a `ValueError` is raised
      to prevent ambiguity in the annealing process.
    - The resulting alignment is calculated with stagger (`ovhg`), ensuring proper
      alignment of the two DNA strands.

    Examples
    --------

    >>> from pydna.dseq import pair
    >>> watson = "ACGTACGT"
    >>> crick =  "TACGTACG"
    >>> pair(watson, crick)
    b'ECGTACGTF'


    """
    watson = watson if isinstance(watson, str) else watson.decode("ascii")
    crick = crick if isinstance(crick, str) else crick.decode("ascii")

    if ovhg is None:

        olaps = _common_sub_strings(
            str(watson).lower(),
            str(_rc(crick).lower()),
            int(_math.log(len(watson)) / _math.log(4)),
        )

        if len(olaps) == 0:
            raise ValueError("Could not anneal the two strands.")

        # We extract the positions and length of the first (longest) overlap, since
        # common_sub_strings sorts the overlaps by length.
        pos_watson, pos_crick, longest_olap_length = olaps[0]

        # We see if there is another overlap of the same length
        if any(olap[2] >= longest_olap_length for olap in olaps[1:]):
            raise ValueError("More than one way of annealing the" " strands.")

        # ovhg is the stagger between the two strings
        ovhg = pos_crick - pos_watson

    # For fragments with a 5' protruding end at the left side

    watson = watson.rjust(len(watson) + ovhg)
    crick = crick[::-1].rjust(len(crick) - ovhg)

    watson = watson.ljust(len(crick))
    crick = crick.ljust(len(watson))

    assert len(watson) == len(crick), "Sense and antisense are not of equal length"

    return "".join(bp_dict_str[pair] for pair in zip(watson, crick)).encode("ascii")


def helper(dseqs):
    print("assert XXX == (" + ", ".join(f'Dseq("{d._data.decode()}")' for d in dseqs) + ")")


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
