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


import copy as _copy
import itertools as _itertools
import re as _re
import sys as _sys
import math as _math

from pydna.seq import Seq as _Seq
from Bio.Seq import _translate_str, _SeqAbstractBaseClass

from pydna._pretty import pretty_str as _pretty_str
from seguid import ldseguid as _ldseguid
from seguid import cdseguid as _cdseguid

from pydna.utils import rc as _rc
from pydna.utils import flatten as _flatten
from pydna.utils import cuts_overlap as _cuts_overlap

from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
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

if TYPE_CHECKING:
    from Bio.Restriction import AbstractCut as _AbstractCut


# To represent any subclass of Dseq
DseqType = _TypeVar("DseqType", bound="Dseq")
EnzymesType = _TypeVar("EnzymesType", _RestrictionBatch, _Iterable["_AbstractCut"], "_AbstractCut")
CutSiteType = _Tuple[_Tuple[int, int], _Union["_AbstractCut", None]]


class Dseq(_Seq):
    """Dseq holds information for a double stranded DNA fragment.

    Dseq also holds information describing the topology of
    the DNA fragment (linear or circular).

    Parameters
    ----------
    watson : str
        a string representing the watson (sense) DNA strand.

    crick : str, optional
        a string representing the crick (antisense) DNA strand.

    ovhg : int, optional
        A positive or negative number to describe the stagger between the
        watson and crick strands.
        see below for a detailed explanation.

    linear : bool, optional
        True indicates that sequence is linear, False that it is circular.

    circular : bool, optional
        True indicates that sequence is circular, False that it is linear.


    Examples
    --------
    Dseq is a subclass of the Biopython Seq object. It stores two
    strings representing the watson (sense) and crick(antisense) strands.
    two properties called linear and circular, and a numeric value ovhg
    (overhang) describing the stagger for the watson and crick strand
    in the 5' end of the fragment.

    The most common usage is probably to create a Dseq object as a
    part of a Dseqrecord object (see :class:`pydna.dseqrecord.Dseqrecord`).

    There are three ways of creating a Dseq object directly listed below, but you can also
    use the function Dseq.from_full_sequence_and_overhangs() to create a Dseq:

    Only one argument (string):

    >>> from pydna.dseq import Dseq
    >>> Dseq("aaa")
    Dseq(-3)
    aaa
    ttt

    The given string will be interpreted as the watson strand of a
    blunt, linear double stranded sequence object. The crick strand
    is created automatically from the watson strand.

    Two arguments (string, string):

    >>> from pydna.dseq import Dseq
    >>> Dseq("gggaaat","ttt")
    Dseq(-7)
    gggaaat
       ttt

    If both watson and crick are given, but not ovhg an attempt
    will be made to find the best annealing between the strands.
    There are limitations to this. For long fragments it is quite
    slow. The length of the annealing sequences have to be at least
    half the length of the shortest of the strands.

    Three arguments (string, string, ovhg=int):

    The ovhg parameter is an integer describing the length of the
    crick strand overhang in the 5' end of the molecule.

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

    >>> Dseq(watson="agt", crick="actta", ovhg=-2)
    Dseq(-7)
    agt
      attca
    >>> Dseq(watson="agt",crick="actta",ovhg=-1)
    Dseq(-6)
    agt
     attca
    >>> Dseq(watson="agt",crick="actta",ovhg=0)
    Dseq(-5)
    agt
    attca
    >>> Dseq(watson="agt",crick="actta",ovhg=1)
    Dseq(-5)
     agt
    attca
    >>> Dseq(watson="agt",crick="actta",ovhg=2)
    Dseq(-5)
      agt
    attca

    If the ovhg parameter is specified a crick strand also
    needs to be supplied, otherwise an exception is raised.

    >>> Dseq(watson="agt", ovhg=2)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/dist-packages/pydna_/dsdna.py", line 169, in __init__
        else:
    ValueError: ovhg defined without crick strand!


    The shape of the fragment is set by circular = True, False

    Note that both ends of the DNA fragment has to be compatible to set
    circular = True.


    >>> Dseq("aaa","ttt")
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=0)
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=1)
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

    trunc = 30

    def __init__(
        self,
        watson: _Union[str, bytes],
        crick: _Union[str, bytes, None] = None,
        ovhg=None,
        circular=False,
        pos=0,
    ):
        if isinstance(watson, bytes):
            watson = watson.decode("ASCII")
        if isinstance(crick, bytes):
            crick = crick.decode("ASCII")

        if crick is None:
            if ovhg is not None:
                raise ValueError("ovhg defined without crick strand!")
            crick = _rc(watson)
            ovhg = 0
            self._data = bytes(watson, encoding="ASCII")

        else:  # crick strand given
            if ovhg is None:  # ovhg not given
                olaps = _common_sub_strings(
                    str(watson).lower(),
                    str(_rc(crick).lower()),
                    int(_math.log(len(watson)) / _math.log(4)),
                )
                if len(olaps) == 0:
                    raise ValueError("Could not anneal the two strands." " Please provide ovhg value")

                # We extract the positions and length of the first (longest) overlap, since
                # common_sub_strings sorts the overlaps by length.
                pos_watson, pos_crick, longest_olap_length = olaps[0]

                # We see if there is another overlap of the same length
                if any(olap[2] >= longest_olap_length for olap in olaps[1:]):
                    raise ValueError("More than one way of annealing the" " strands. Please provide ovhg value")

                ovhg = pos_crick - pos_watson

                sns = (ovhg * " ") + _pretty_str(watson)
                asn = (-ovhg * " ") + _pretty_str(_rc(crick))

                self._data = bytes(
                    "".join([a.strip() or b.strip() for a, b in _itertools.zip_longest(sns, asn, fillvalue=" ")]),
                    encoding="ASCII",
                )

            else:  # ovhg given
                if ovhg == 0:
                    if len(watson) >= len(crick):
                        self._data = bytes(watson, encoding="ASCII")
                    else:
                        self._data = bytes(
                            watson + _rc(crick[: len(crick) - len(watson)]),
                            encoding="ASCII",
                        )
                elif ovhg > 0:
                    if ovhg + len(watson) > len(crick):
                        self._data = bytes(_rc(crick[-ovhg:]) + watson, encoding="ASCII")
                    else:
                        self._data = bytes(
                            _rc(crick[-ovhg:]) + watson + _rc(crick[: len(crick) - ovhg - len(watson)]),
                            encoding="ASCII",
                        )
                else:  # ovhg < 0
                    if -ovhg + len(crick) > len(watson):
                        self._data = bytes(
                            watson + _rc(crick[: -ovhg + len(crick) - len(watson)]),
                            encoding="ASCII",
                        )
                    else:
                        self._data = bytes(watson, encoding="ASCII")

        self.circular = circular
        self.watson = _pretty_str(watson)
        self.crick = _pretty_str(crick)
        self.length = len(self._data)
        self.ovhg = ovhg
        self.pos = pos

    @classmethod
    def quick(
        cls,
        watson: str,
        crick: str,
        ovhg=0,
        circular=False,
        pos=0,
    ):
        obj = cls.__new__(cls)  # Does not call __init__
        obj.watson = _pretty_str(watson)
        obj.crick = _pretty_str(crick)
        obj.ovhg = ovhg
        obj.circular = circular
        obj.length = max(len(watson) + max(0, ovhg), len(crick) + max(0, -ovhg))
        obj.pos = pos
        wb = bytes(watson, encoding="ASCII")
        cb = bytes(crick, encoding="ASCII")
        obj._data = _rc(cb[-max(0, ovhg) or len(cb) :]) + wb + _rc(cb[: max(0, len(cb) - ovhg - len(wb))])
        return obj

    @classmethod
    def from_string(
        cls,
        dna: str,
        *args,
        # linear=True,
        circular=False,
        **kwargs,
    ):
        obj = cls.__new__(cls)  # Does not call __init__
        obj.watson = _pretty_str(dna)
        obj.crick = _pretty_str(_rc(dna))
        obj.ovhg = 0
        obj.circular = circular
        # obj._linear = linear
        obj.length = len(dna)
        obj.pos = 0
        obj._data = bytes(dna, encoding="ASCII")
        return obj

    @classmethod
    def from_representation(cls, dsdna: str, *args, **kwargs):
        obj = cls.__new__(cls)  # Does not call __init__
        w, c, *r = [ln for ln in dsdna.splitlines() if ln]
        ovhg = obj.ovhg = len(w) - len(w.lstrip()) - (len(c) - len(c.lstrip()))
        watson = obj.watson = _pretty_str(w.strip())
        crick = obj.crick = _pretty_str(c.strip()[::-1])
        obj.circular = False
        # obj._linear = True
        obj.length = max(len(watson) + max(0, ovhg), len(crick) + max(0, -ovhg))
        obj.pos = 0
        wb = bytes(watson, encoding="ASCII")
        cb = bytes(crick, encoding="ASCII")
        obj._data = _rc(cb[-max(0, ovhg) or len(cb) :]) + wb + _rc(cb[: max(0, len(cb) - ovhg - len(wb))])
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

    # @property
    # def ovhg(self):
    #     """The ovhg property. This cannot be set directly, but is a
    #     consequence of how the watson and crick strands anneal to
    #     each other"""
    #     return self._ovhg

    # @property
    # def linear(self):
    #     """The linear property can not be set directly.
    #     Use an empty slice [:] to create a linear object."""
    #     return self._linear

    # @property
    # def circular(self):
    #     """The circular property can not be set directly.
    #     Use :meth:`looped` to create a circular Dseq object"""
    #     return self._circular

    def mw(self) -> float:
        """This method returns the molecular weight of the DNA molecule
        in g/mol. The following formula is used::

               MW = (A x 313.2) + (T x 304.2) +
                    (C x 289.2) + (G x 329.2) +
                    (N x 308.9) + 79.0
        """
        nts = (self.watson + self.crick).lower()

        return (
            313.2 * nts.count("a")
            + 304.2 * nts.count("t")
            + 289.2 * nts.count("c")
            + 329.2 * nts.count("g")
            + 308.9 * nts.count("n")
            + 79.0
        )

    def upper(self: DseqType) -> DseqType:
        """Return an upper case copy of the sequence.

        >>> from pydna.dseq import Dseq
        >>> my_seq = Dseq("aAa")
        >>> my_seq
        Dseq(-3)
        aAa
        tTt
        >>> my_seq.upper()
        Dseq(-3)
        AAA
        TTT

        Returns
        -------
        Dseq
            Dseq object in uppercase

        See also
        --------
        pydna.dseq.Dseq.lower

        """
        return self.quick(
            self.watson.upper(),
            self.crick.upper(),
            ovhg=self.ovhg,
            # linear=self.linear,
            circular=self.circular,
            pos=self.pos,
        )

    def lower(self: DseqType) -> DseqType:
        """Return a lower case copy of the sequence.

        >>> from pydna.dseq import Dseq
        >>> my_seq = Dseq("aAa")
        >>> my_seq
        Dseq(-3)
        aAa
        tTt
        >>> my_seq.lower()
        Dseq(-3)
        aaa
        ttt

        Returns
        -------
        Dseq
            Dseq object in lowercase

        See also
        --------
        pydna.dseq.Dseq.upper
        """
        return self.quick(
            self.watson.lower(),
            self.crick.lower(),
            ovhg=self.ovhg,
            # linear=self.linear,
            circular=self.circular,
            pos=self.pos,
        )

    def find(self, sub: _Union[_SeqAbstractBaseClass, str, bytes], start=0, end=_sys.maxsize) -> int:
        """This method behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Returns -1 if the subsequence is NOT found.

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
        >>> seq = Dseq("atcgactgacgtgtt")
        >>> seq
        Dseq(-15)
        atcgactgacgtgtt
        tagctgactgcacaa
        >>> seq.find("gac")
        3
        >>> seq = Dseq(watson="agt",crick="actta",ovhg=-2)
        >>> seq
        Dseq(-7)
        agt
          attca
        >>> seq.find("taa")
        2
        """

        if not self.circular:
            return _Seq.find(self, sub, start, end)

        return (_pretty_str(self) + _pretty_str(self)).find(sub, start, end)

    def __getitem__(self, sl: slice) -> "Dseq":
        """Returns a subsequence. This method is used by the slice notation"""

        if not self.circular:
            x = len(self.crick) - self.ovhg - len(self.watson)

            sns = (self.ovhg * " " + self.watson + x * " ")[sl]
            asn = (-self.ovhg * " " + self.crick[::-1] + -x * " ")[sl]

            ovhg = max((len(sns) - len(sns.lstrip()), -len(asn) + len(asn.lstrip())), key=abs)

            return Dseq(
                sns.strip(),
                asn[::-1].strip(),
                ovhg=ovhg,
                # linear=True
            )
        else:
            sl = slice(sl.start or 0, sl.stop or len(self), sl.step)
            if sl.start > len(self) or sl.stop > len(self):
                return Dseq("")
            if sl.start < sl.stop:
                return Dseq(
                    self.watson[sl],
                    self.crick[::-1][sl][::-1],
                    ovhg=0,
                    # linear=True
                )
            else:
                try:
                    stp = abs(sl.step)
                except TypeError:
                    stp = 1
                start = sl.start
                stop = sl.stop

                w = self.watson[(start or len(self)) :: stp] + self.watson[: (stop or 0) : stp]
                c = self.crick[len(self) - stop :: stp] + self.crick[: len(self) - start : stp]

                return Dseq(w, c, ovhg=0)  # , linear=True)

    def __eq__(self, other: DseqType) -> bool:
        """Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties. This comparison is case
        insensitive.

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

    def __repr__(self):
        """Returns a representation of the sequence, truncated if
        longer than 30 bp"""

        if len(self) > Dseq.trunc:
            if self.ovhg > 0:
                d = self.crick[-self.ovhg :][::-1]
                hej = len(d)
                if len(d) > 10:
                    d = "{}..{}".format(d[:4], d[-4:])
                a = len(d) * " "

            elif self.ovhg < 0:
                a = self.watson[: max(0, -self.ovhg)]
                hej = len(a)
                if len(a) > 10:
                    a = "{}..{}".format(a[:4], a[-4:])
                d = len(a) * " "
            else:
                a = ""
                d = ""
                hej = 0

            x = self.ovhg + len(self.watson) - len(self.crick)

            if x > 0:
                c = self.watson[len(self.crick) - self.ovhg :]
                y = len(c)
                if len(c) > 10:
                    c = "{}..{}".format(c[:4], c[-4:])
                f = len(c) * " "
            elif x < 0:
                f = self.crick[:-x][::-1]
                y = len(f)
                if len(f) > 10:
                    f = "{}..{}".format(f[:4], f[-4:])
                c = len(f) * " "
            else:
                c = ""
                f = ""
                y = 0

            L = len(self) - hej - y
            x1 = -min(0, self.ovhg)
            x2 = x1 + L
            x3 = -min(0, x)
            x4 = x3 + L

            b = self.watson[x1:x2]
            e = self.crick[x3:x4][::-1]

            if len(b) > 10:
                b = "{}..{}".format(b[:4], b[-4:])
                e = "{}..{}".format(e[:4], e[-4:])

            return _pretty_str("{klass}({top}{size})\n" "{a}{b}{c}\n" "{d}{e}{f}").format(
                klass=self.__class__.__name__,
                top={False: "-", True: "o"}[self.circular],
                size=len(self),
                a=a,
                b=b,
                c=c,
                d=d,
                e=e,
                f=f,
            )

        else:
            return _pretty_str(
                "{}({}{})\n{}\n{}".format(
                    self.__class__.__name__,
                    {False: "-", True: "o"}[self.circular],
                    len(self),
                    self.ovhg * " " + self.watson,
                    -self.ovhg * " " + self.crick[::-1],
                )
            )

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
        return Dseq.quick(
            self.crick,
            self.watson,
            ovhg=len(self.watson) - len(self.crick) + self.ovhg,
            circular=self.circular,
        )

    rc = reverse_complement  # alias for reverse_complement

    def shifted(self: DseqType, shift: int) -> DseqType:
        """Shifted version of a circular Dseq object."""
        if not self.circular:
            raise TypeError("DNA is not circular.")
        shift = shift % len(self)
        if not shift:
            return _copy.deepcopy(self)
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
        >>> a.T4("t")
        Dseq(-8)
        catcgat
         tagctag
        >>> a.T4("t").looped()
        Dseq(o7)
        catcgat
        gtagcta
        >>> a.T4("a")
        Dseq(-8)
        catcga
          agctag
        >>> a.T4("a").looped()
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
        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()
        if type5 == type3 and str(sticky5) == str(_rc(sticky3)):
            nseq = self.__class__.quick(
                self.watson,
                self.crick[-self.ovhg :] + self.crick[: -self.ovhg],
                ovhg=0,
                # linear=False,
                circular=True,
            )
            # assert len(nseq.crick) == len(nseq.watson)
            return nseq
        else:
            raise TypeError("DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!")

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
        return selfcopy  # self.__class__(self.watson, linear=True)

    def five_prime_end(self) -> _Tuple[str, str]:
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.five_prime_end()
        ('blunt', '')
        >>> a=Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.five_prime_end()
        ("3'", 't')
        >>> a=Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.five_prime_end()
        ("5'", 'a')
        >>>

        See also
        --------
        pydna.dseq.Dseq.three_prime_end

        """
        if self.watson and not self.crick:
            return "5'", self.watson.lower()
        if not self.watson and self.crick:
            return "3'", self.crick.lower()
        if self.ovhg < 0:
            sticky = self.watson[: -self.ovhg].lower()
            type_ = "5'"
        elif self.ovhg > 0:
            sticky = self.crick[-self.ovhg :].lower()
            type_ = "3'"
        else:
            sticky = ""
            type_ = "blunt"
        return type_, sticky

    def three_prime_end(self) -> _Tuple[str, str]:
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', '')
        >>> a=Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.three_prime_end()
        ("3'", 'a')
        >>> a=Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.three_prime_end()
        ("5'", 't')
        >>>

        See also
        --------
        pydna.dseq.Dseq.five_prime_end

        """

        ovhg = len(self.watson) - len(self.crick) + self.ovhg

        if ovhg < 0:
            sticky = self.crick[:-ovhg].lower()
            type_ = "5'"
        elif ovhg > 0:
            sticky = self.watson[-ovhg:].lower()
            type_ = "3'"
        else:
            sticky = ""
            type_ = "blunt"
        return type_, sticky

    def watson_ovhg(self) -> int:
        """Returns the overhang of the watson strand at the three prime."""
        return len(self.watson) - len(self.crick) + self.ovhg

    def __add__(self: DseqType, other: DseqType) -> DseqType:
        """Simulates ligation between two DNA fragments.

        Add other Dseq object at the end of the sequence.
        Type error is raised if any of the points below are fulfilled:

        * one or more objects are circular
        * if three prime sticky end of self is not the same type
          (5' or 3') as the sticky end of other
        * three prime sticky end of self complementary with five
          prime sticky end of other.

        Phosphorylation and dephosphorylation is not considered.

        DNA is allways presumed to have the necessary 5' phospate
        group necessary for ligation.

        """
        # test for circular DNA
        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        self_type, self_tail = self.three_prime_end()
        other_type, other_tail = other.five_prime_end()

        if self_type == other_type and str(self_tail) == str(_rc(other_tail)):
            answer = Dseq.quick(self.watson + other.watson, other.crick + self.crick, self.ovhg)
        elif not self:
            answer = _copy.deepcopy(other)
        elif not other:
            answer = _copy.deepcopy(self)
        else:
            raise TypeError("sticky ends not compatible!")
        return answer

    def __mul__(self: DseqType, number: int) -> DseqType:
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(number)))
        if number <= 0:
            return self.__class__("")
        new = _copy.deepcopy(self)
        for i in range(number - 1):
            new += self
        return new

    def _fill_in_five_prime(self: DseqType, nucleotides: str) -> str:
        stuffer = ""
        type, se = self.five_prime_end()
        if type == "5'":
            for n in _rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.crick + stuffer, self.ovhg + len(stuffer)

    def _fill_in_three_prime(self: DseqType, nucleotides: str) -> str:
        stuffer = ""
        type, se = self.three_prime_end()
        if type == "5'":
            for n in _rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.watson + stuffer

    def fill_in(self, nucleotides: _Union[None, str] = None) -> "Dseq":
        """Fill in of five prime protruding end with a DNA polymerase
        that has only DNA polymerase activity (such as exo-klenow [#]_)
        and any combination of A, G, C or T. Default are all four
        nucleotides together.

        Parameters
        ----------

        nucleotides : str

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.fill_in()
        Dseq(-3)
        aaa
        ttt
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
        >>>

        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Klenow_fragment#The_exo-_Klenow_fragment

        """
        if nucleotides is None:
            nucleotides = "GATCRYWSMKHBVDN"

        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        crick, ovhg = self._fill_in_five_prime(nucleotides)
        watson = self._fill_in_three_prime(nucleotides)
        return Dseq(watson, crick, ovhg)

    def transcribe(self) -> _Seq:
        return _Seq(self.watson).transcribe()

    def translate(self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-") -> _Seq:
        return _Seq(_translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap))

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
        return Dseq(self.watson[max(0, -self.ovhg) : min(len(self.watson), len(self.crick) - self.ovhg)])

    def T4(self, nucleotides=None) -> "Dseq":
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
        >>> a.T4()
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

        if not nucleotides:
            nucleotides = "GATCRYWSMKHBVDN"
        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        type, se = self.five_prime_end()
        if type == "5'":
            crick, ovhg = self._fill_in_five_prime(nucleotides)
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
            watson = self._fill_in_three_prime(nucleotides)
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
            cs = _cdseguid(self.watson.upper(), self.crick.upper(), alphabet="{DNA-extended}")
        else:
            """docstring."""
            w = f"{self.ovhg * '-'}{self.watson}{'-' * (-self.ovhg + len(self.crick) - len(self.watson))}".upper()
            c = f"{'-' * (self.ovhg + len(self.watson) - len(self.crick))}{self.crick}{-self.ovhg * '-'}".upper()
            cs = _ldseguid(w, c, alphabet="{DNA-extended}")
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
        return self.ovhg == 0 and len(self.watson) == len(self.crick) and not self.circular

    def cas9(self, RNA: str) -> _Tuple[slice, ...]:
        """docstring."""
        bRNA = bytes(RNA, "ASCII")
        slices = []
        cuts = [0]
        for m in _re.finditer(bRNA, self._data):
            cuts.append(m.start() + 17)
        cuts.append(self.length)
        slices = tuple(slice(x, y, 1) for x, y in zip(cuts, cuts[1:]))
        return slices

    def terminal_transferase(self, nucleotides="a") -> "Dseq":
        """docstring."""
        ovhg = self.ovhg
        if self.ovhg >= 0:
            ovhg += len(nucleotides)
        return Dseq(self.watson + nucleotides, self.crick + nucleotides, ovhg)

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
        cutsite_pairs = self.get_cutsite_pairs(cutsites)
        return tuple(self.apply_cut(*cs) for cs in cutsite_pairs)

    def cutsite_is_valid(self, cutsite: CutSiteType) -> bool:
        """Returns False if:
        - Cut positions fall outside the sequence (could be moved to Biopython)
        - Overhang is not double stranded
        - Recognition site is not double stranded or is outside the sequence
        - For enzymes that cut twice, it checks that at least one possibility is valid
        """

        assert cutsite is not None, "cutsite is None"

        enz = cutsite[1]
        watson, crick, ovhg = self.get_cut_parameters(cutsite, True)

        # The overhang is double stranded
        overhang_dseq = self[watson:crick] if ovhg < 0 else self[crick:watson]
        if overhang_dseq.ovhg != 0 or overhang_dseq.watson_ovhg() != 0:
            return False

        # The recognition site is double stranded and within the sequence
        start_of_recognition_site = watson - enz.fst5
        if start_of_recognition_site < 0:
            start_of_recognition_site += len(self)
        end_of_recognition_site = start_of_recognition_site + enz.size
        if self.circular:
            end_of_recognition_site %= len(self)
        recognition_site = self[start_of_recognition_site:end_of_recognition_site]
        if len(recognition_site) == 0 or recognition_site.ovhg != 0 or recognition_site.watson_ovhg() != 0:
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
                recognition_site = self[start_of_recognition_site:end_of_recognition_site]

                if len(recognition_site) == 0 or recognition_site.ovhg != 0 or recognition_site.watson_ovhg() != 0:
                    return False

        return True

    def get_cutsites(self: DseqType, *enzymes: EnzymesType) -> _List[CutSiteType]:
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
        out = list()
        for e in enzymes:
            # Positions of the cut on the watson strand. They are 1-based, so we subtract
            # 1 to get 0-based positions
            cuts_watson = [c - 1 for c in e.search(self, linear=(not self.circular))]

            out += [((w, e.ovhg), e) for w in cuts_watson]

        return sorted([cutsite for cutsite in out if self.cutsite_is_valid(cutsite)])

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

    def get_cut_parameters(self, cut: _Union[CutSiteType, None], is_left: bool) -> _Tuple[int, int, int]:
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
        return *self.right_end_position(), self.watson_ovhg()

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
        if _cuts_overlap(left_cut, right_cut, len(self)):
            raise ValueError("Cuts by {} {} overlap.".format(left_cut[1], right_cut[1]))

        left_watson, left_crick, ovhg_left = self.get_cut_parameters(left_cut, True)
        right_watson, right_crick, _ = self.get_cut_parameters(right_cut, False)
        return Dseq(
            str(self[left_watson:right_watson]),
            # The line below could be easier to understand as _rc(str(self[left_crick:right_crick])), but it does not preserve the case
            str(self.reverse_complement()[len(self) - right_crick : len(self) - left_crick]),
            ovhg=ovhg_left,
        )

    def get_cutsite_pairs(
        self, cutsites: _List[CutSiteType]
    ) -> _List[_Tuple[_Union[None, CutSiteType], _Union[None, CutSiteType]]]:
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


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
