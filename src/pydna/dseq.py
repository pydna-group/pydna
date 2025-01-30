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


class oldDseq(_Seq):
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
        full_sequence_rev = str(oldDseq(full_sequence).reverse_complement())
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

        return oldDseq(watson, crick=crick, ovhg=crick_ovhg)

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
        >>> my_seq = oldDseq("aAa")
        >>> my_seq
        oldDseq(-3)
        aAa
        tTt
        >>> my_seq.upper()
        oldDseq(-3)
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
        >>> my_seq = oldDseq("aAa")
        >>> my_seq
        oldDseq(-3)
        aAa
        tTt
        >>> my_seq.lower()
        oldDseq(-3)
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
        >>> seq = oldDseq("atcgactgacgtgtt")
        >>> seq
        oldDseq(-15)
        atcgactgacgtgtt
        tagctgactgcacaa
        >>> seq.find("gac")
        3
        >>> seq = oldDseq(watson="agt",crick="actta",ovhg=-2)
        >>> seq
        oldDseq(-7)
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

            return oldDseq(
                sns.strip(),
                asn[::-1].strip(),
                ovhg=ovhg,
                # linear=True
            )
        else:
            sl = slice(sl.start or 0, sl.stop or len(self), sl.step)
            if sl.start > len(self) or sl.stop > len(self):
                return oldDseq("")
            if sl.start < sl.stop:
                return oldDseq(
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

                return oldDseq(w, c, ovhg=0)  # , linear=True)

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
        >>> a=oldDseq("catcgatc")
        >>> a
        oldDseq(-8)
        catcgatc
        gtagctag
        >>> b=a.reverse_complement()
        >>> b
        oldDseq(-8)
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
        >>> a=oldDseq("catcgatc")
        >>> a
        oldDseq(-8)
        catcgatc
        gtagctag
        >>> a.looped()
        oldDseq(o8)
        catcgatc
        gtagctag
        >>> a.T4("t")
        oldDseq(-8)
        catcgat
         tagctag
        >>> a.T4("t").looped()
        oldDseq(o7)
        catcgat
        gtagcta
        >>> a.T4("a")
        oldDseq(-8)
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
        >>> a=oldDseq("catcgatc", circular=True)
        >>> a
        oldDseq(o8)
        catcgatc
        gtagctag
        >>> a[:]
        oldDseq(-8)
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
        >>> a=oldDseq("aaa", "ttt")
        >>> a
        oldDseq(-3)
        aaa
        ttt
        >>> a.five_prime_end()
        ('blunt', '')
        >>> a=oldDseq("aaa", "ttt", ovhg=1)
        >>> a
        oldDseq(-4)
         aaa
        ttt
        >>> a.five_prime_end()
        ("3'", 't')
        >>> a=oldDseq("aaa", "ttt", ovhg=-1)
        >>> a
        oldDseq(-4)
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
        >>> a=oldDseq("aaa", "ttt")
        >>> a
        oldDseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', '')
        >>> a=oldDseq("aaa", "ttt", ovhg=1)
        >>> a
        oldDseq(-4)
         aaa
        ttt
        >>> a.three_prime_end()
        ("3'", 'a')
        >>> a=oldDseq("aaa", "ttt", ovhg=-1)
        >>> a
        oldDseq(-4)
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
        >>> a=oldDseq("aaa", "ttt")
        >>> a
        oldDseq(-3)
        aaa
        ttt
        >>> a.fill_in()
        oldDseq(-3)
        aaa
        ttt
        >>> b=oldDseq("caaa", "cttt")
        >>> b
        oldDseq(-5)
        caaa
         tttc
        >>> b.fill_in()
        oldDseq(-5)
        caaag
        gtttc
        >>> b.fill_in("g")
        oldDseq(-5)
        caaag
        gtttc
        >>> b.fill_in("tac")
        oldDseq(-5)
        caaa
         tttc
        >>> c=oldDseq("aaac", "tttg")
        >>> c
        oldDseq(-5)
         aaac
        gttt
        >>> c.fill_in()
        oldDseq(-5)
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
        return oldDseq(watson, crick, ovhg)

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
         >>> b=oldDseq("caaa", "cttt")
         >>> b
         oldDseq(-5)
         caaa
          tttc
         >>> b.mung()
         oldDseq(-3)
         aaa
         ttt
         >>> c=oldDseq("aaac", "tttg")
         >>> c
         oldDseq(-5)
          aaac
         gttt
         >>> c.mung()
         oldDseq(-3)
         aaa
         ttt



        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Mung_bean_nuclease


        """
        return oldDseq(self.watson[max(0, -self.ovhg) : min(len(self.watson), len(self.crick) - self.ovhg)])

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
        >>> a=oldDseq("gatcgatc")
        >>> a
        oldDseq(-8)
        gatcgatc
        ctagctag
        >>> a.T4()
        oldDseq(-8)
        gatcgatc
        ctagctag
        >>> a.T4("t")
        oldDseq(-8)
        gatcgat
         tagctag
        >>> a.T4("a")
        oldDseq(-8)
        gatcga
          agctag
        >>> a.T4("g")
        oldDseq(-8)
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
        return oldDseq(watson, crick, ovhg)

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
        >>> a=oldDseq("gat")
        >>> a
        oldDseq(-3)
        gat
        cta
        >>> a.isblunt()
        True
        >>> a=oldDseq("gat", "atcg")
        >>> a
        oldDseq(-4)
         gat
        gcta
        >>> a.isblunt()
        False
        >>> a=oldDseq("gat", "gatc")
        >>> a
        oldDseq(-4)
        gat
        ctag
        >>> a.isblunt()
        False
        >>> a=oldDseq("gat", circular=True)
        >>> a
        oldDseq(o3)
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
        return oldDseq(self.watson + nucleotides, self.crick + nucleotides, ovhg)

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
        >>> seq=oldDseq("ggatccnnngaattc")
        >>> seq
        oldDseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>> from Bio.Restriction import BamHI,EcoRI
        >>> type(seq.cut(BamHI))
        <class 'tuple'>
        >>> for frag in seq.cut(BamHI): print(repr(frag))
        oldDseq(-5)
        g
        cctag
        oldDseq(-14)
        gatccnnngaattc
            gnnncttaag
        >>> seq.cut(EcoRI, BamHI) ==  seq.cut(BamHI, EcoRI)
        True
        >>> a,b,c = seq.cut(EcoRI, BamHI)
        >>> a+b+c
        oldDseq(-15)
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
        >>> seq = oldDseq('AAGAATTCAAGAATTC')
        >>> seq.get_cutsites(EcoRI)
        [((3, -4), EcoRI), ((11, -4), EcoRI)]

        `cut_watson` is defined with respect to the "full sequence", not the
        watson strand:

        >>> dseq = Dseq.from_full_sequence_and_overhangs('aaGAATTCaa', 1, 0)
        >>> dseq
        oldDseq(-10)
         aGAATTCaa
        ttCTTAAGtt
        >>> dseq.get_cutsites([EcoRI])
        [((3, -4), EcoRI)]

        Cuts are only returned if the recognition site and overhang are on the double-strand
        part of the sequence.

        >>> oldDseq('GAATTC').get_cutsites([EcoRI])
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
        >>> dseq = oldDseq('aaGAATTCaaGAATTCaa')
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((3, -4), EcoRI), ((11, -4), EcoRI)]
        >>> p1, p2, p3 = dseq.get_cutsite_pairs(cutsites)
        >>> p1
        (None, ((3, -4), EcoRI))
        >>> dseq.apply_cut(*p1)
        oldDseq(-7)
        aaG
        ttCTTAA
        >>> p2
        (((3, -4), EcoRI), ((11, -4), EcoRI))
        >>> dseq.apply_cut(*p2)
        oldDseq(-12)
        AATTCaaG
            GttCTTAA
        >>> p3
        (((11, -4), EcoRI), None)
        >>> dseq.apply_cut(*p3)
        oldDseq(-7)
        AATTCaa
            Gtt

        >>> dseq = oldDseq('TTCaaGAA', circular=True)
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((6, -4), EcoRI)]
        >>> pair = dseq.get_cutsite_pairs(cutsites)[0]
        >>> pair
        (((6, -4), EcoRI), ((6, -4), EcoRI))
        >>> dseq.apply_cut(*pair)
        oldDseq(-12)
        AATTCaaG
            GttCTTAA

        """
        if _cuts_overlap(left_cut, right_cut, len(self)):
            raise ValueError("Cuts by {} {} overlap.".format(left_cut[1], right_cut[1]))

        left_watson, left_crick, ovhg_left = self.get_cut_parameters(left_cut, True)
        right_watson, right_crick, _ = self.get_cut_parameters(right_cut, False)
        print(left_watson, left_crick, ovhg_left)
        print(right_watson, right_crick, _)
        return oldDseq(
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
        >>> dseq = oldDseq('aaGAATTCaaGAATTCaa')
        >>> cutsites = dseq.get_cutsites([EcoRI])
        >>> cutsites
        [((3, -4), EcoRI), ((11, -4), EcoRI)]
        >>> dseq.get_cutsite_pairs(cutsites)
        [(None, ((3, -4), EcoRI)), (((3, -4), EcoRI), ((11, -4), EcoRI)), (((11, -4), EcoRI), None)]

        >>> dseq = oldDseq('TTCaaGAA', circular=True)
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
from Bio.SeqFeature import SimpleLocation as _lc

from seguid import ldseguid as _ldseguid
from seguid import cdseguid as _cdseguid

from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import shift_location as _sl
from pydna.utils import rc as _rc
from pydna.utils import flatten as _flatten

from pydna.utils import cuts_overlap2 as _cuts_overlap2
from pydna.utils import complement as _complement
from pydna.utils import bp_dict
from pydna.utils import bp_dict_str
from pydna.utils import to_watson_table
from pydna.utils import to_crick_table
from pydna.utils import to_5tail_table
from pydna.utils import to_3tail_table
from pydna.utils import to_full_sequence
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
        cutsite_pairs = self.get_cutsite_pairs(cutsites)
        return tuple(self.apply_cut(*cutsite_pair) for cutsite_pair in cutsite_pairs)

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
            # if _cuts_overlap(left_cut, right_cut, len(self)):
            #     raise ValueError("Cuts by {} {} overlap.".format(left_cut[1], right_cut[1]))
            # left_watson, left_crick, ovhg_left = self.get_cut_parameters(left_cut, True)
            # right_watson, right_crick, _ = self.get_cut_parameters(right_cut, False)

            (left_begin, left_overhang), *rest = left_cut
            (right_begin, right_overhang), *rest = right_cut

            table1 = to_5tail_table if left_overhang>0 else to_3tail_table
            table2 = to_5tail_table if right_overhang<0 else to_3tail_table

            left_stck_begin = min(left_begin, left_begin - left_overhang)
            left_stck_end = left_stck_begin + abs(left_overhang)

            right_stck_begin = min(right_begin, right_begin - right_overhang)
            right_stck_end = right_stck_begin + abs(right_overhang)

            return Dseq(
                self._data[left_stck_begin:left_stck_end].translate(table1)
                + self._data[left_stck_end:right_stck_begin]
                + self._data[right_stck_begin:right_stck_end].translate(table2)
            )

            # return Dseq(
            #     str(self[left_watson:right_watson]),
            #     # The line below could be easier to understand as _rc(str(self[left_crick:right_crick])), but it does not preserve the case
            #     str(self.reverse_complement()[len(self) - right_crick : len(self) - left_crick]),
            #     ovhg=ovhg_left,
            # )

    def get_cutsite_pairs(self, cutsites: _List[CutSiteType]) -> _List[_Tuple[_Union[None, CutSiteType], _Union[None, CutSiteType]]]:
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
        [((3, -4), EcoRI, SimpleLocation(ExactPosition(2), ExactPosition(8))),
         ((11, -4), EcoRI, SimpleLocation(ExactPosition(10), ExactPosition(16)))]

        `cut_watson` is defined with respect to the "full sequence", not the
        watson strand:

        >>> dseq = Dseq.from_full_sequence_and_overhangs('aaGAATTCaa', 1, 0)
        >>> dseq
        Dseq(-10)
         aGAATTCaa
        ttCTTAAGtt
        >>> dseq.get_cutsites([EcoRI])
        [((3, -4), EcoRI, SimpleLocation(ExactPosition(2), ExactPosition(8)))]

        Cuts are only returned if the recognition site and overhang are on the double-strand
        part of the sequence.

        >>> Dseq('GAATTC').get_cutsites([EcoRI])
        [((1, -4), EcoRI, SimpleLocation(ExactPosition(0), ExactPosition(6)))]
        >>> Dseq.from_full_sequence_and_overhangs('GAATTC', -1, 0).get_cutsites([EcoRI])
        []

        """
        if len(enzymes) == 1 and isinstance(enzymes[0], _RestrictionBatch):
            # argument is probably a RestrictionBatch
            enzymes = [e for e in enzymes[0]]

        enzymes = _flatten(enzymes)
        out = list()

        search = str(self).upper()
        ln = len(self)
        # This algorithm for finding restriction sites uses the compsite property of the restriction
        # enzyme instead of the enzyme search method. This is because the Bio.Restriction module
        # coerces the seqquence to a Bio.Restriction.FormattedSeq object and this checks for
        # non-iupac letters in the input. This will not work with pydna since dsiupac is used.
        # https://github.com/MetabolicEngineeringGroupCBMA/MetabolicEngineeringGroupCBMA.github.io/wiki/dsIUPAC
        if self.circular:
            for e in enzymes:
                cuts = []
                csearch = search + search[: len(e) - 1]
                for i in _re.finditer(e.compsite, csearch):
                    on_watson, on_crick, *_ = i.groups() + (None,)
                    if on_watson:
                        cutposition = (i.start() + e.fst5) % ln
                        recsite = _sl(_lc(i.start(), i.start() + e.size), 0, ln)
                        cuts.append(((cutposition, e.ovhg), e, recsite))
                        if e.scd5:
                            cutposition = (i.start() + e.scd5) % ln
                            recsite = _sl(_lc(i.start(), i.start() + e.size), 0, ln)
                            cuts.append(((cutposition, e.ovhg), e, recsite))
                    elif on_crick:
                        cutposition = (i.start() - e.fst3) % ln
                        recsite = _sl(_lc(i.start() - e.size, i.start()), 0, ln)
                        cuts.append(((cutposition, e.ovhg), e, recsite))
                        if e.scd3:
                            cutposition = (i.start() - e.scd3) % ln
                            recsite = _sl(_lc(i.start() - e.size, i.start()), 0, ln)
                            cuts.append(((cutposition, e.ovhg), e, recsite))
                out.extend(cuts)
        else:
            for e in enzymes:
                cuts = []
                for i in _re.finditer(e.compsite, search):
                    on_watson, on_crick, *_ = i.groups() + (None,)

                    if on_watson and i.start() + e.fst5 <= ln and (i.start() + e.fst5 - e.ovhg) <= ln:
                        cutposition = i.start() + e.fst5
                        recsite = _lc(i.start(), i.start() + e.size)
                        cuts.append(((cutposition, e.ovhg), e, recsite))
                    elif on_crick and i.start() - e.fst3 <= ln and (i.start() - e.fst3 + e.ovhg) <= ln:
                        cutposition = i.start() - e.fst3
                        recsite = _lc(i.start() - e.size, i.start())
                        cuts.append(((cutposition, e.ovhg), e, recsite))
                    if e.scd5:
                        if on_watson and i.start() + e.scd5 <= ln and (i.start() + e.scd5 - e.ovhg) <= ln:
                            cutposition = (i.start() + e.scd5)
                            recsite = _lc(i.start(), i.start() + e.size)
                            cuts.append(((cutposition, e.ovhg), e, recsite))
                        elif on_crick and i.start() - e.scd3 <= ln and (i.start() - e.scd3 + e.ovhg) <= ln:
                            cutposition = (i.start() - e.scd3)
                            recsite = _lc(i.start() - e.size, i.start())
                            cuts.append(((cutposition, e.ovhg), e, recsite))
                out.extend(cuts)
        return sorted(out)

    def melt(self: DseqType, length) -> _List[CutSiteType]:
        frags, cutpairs, shift, ln = self._melt_w_params(length)
        return frags

    def _melt_w_params(self: DseqType, length) -> _List[CutSiteType]:

        if length < 1:
            return ()

        frags, meltpairs, shift = [], [], 0

        regex = _re.compile(
            "(^|[PEXIpexiQFZJqfzj])" f"([GATCgatcUuOo]{{1,{length}}})" "([PEXIpexiQFZJqfzj]|$)".encode("ascii")
        )

        def melt(fragment):

            m = regex.search(fragment)

            if not m:
                return None

            before, dsreg, after = m.groups()

            upper = dsreg.translate(to_3tail_table)
            lower = dsreg.translate(to_5tail_table)

            left = cutfrom[: m.start()]
            right = cutfrom[m.end() :]

            if before in b"PEXIpexi" and after in b"QFZJqfzj":
                left = left + before + upper
                right = lower + after + right
            elif before in b"QFZJqfzj" and after in b"PEXIpexi":
                left = left + before + lower
                right = upper + after + right
            elif before in b"PEXIpexi" and after in b"PEXIpexi":
                right = left + before + upper + after + right
                left = lower
            elif before in b"QFZJqfzj" and after in b"QFZJqfzj":
                right = left + before + lower + after + right
                left = upper

            return left, right, 0, m.end(2), m.start(2), len(fragment)

        cutfrom = self._data

        frags = []

        right = b""

        rightstart = 0

        if self.circular:

            left, right, x1, y1, x2, y2 = melt(cutfrom)

            right = cutfrom = right + left

            shift = x2

            meltpairs.append((0, y1))

        while params := melt(cutfrom):

            left, right, x1, y1, x2, y2 = params

            frags.append(left)

            rightstart += x2

            cutfrom = right

            meltpairs.append((x1, y1))

        frags.append(right)

        frags = tuple(Dseq(f) for f in frags if f)

        if not frags:

            meltpairs = ()

        else:

            meltpairs = tuple(meltpairs + [(x2, y2)])

        return frags, meltpairs, shift, len(cutfrom)

    def cutsite_is_valid(self, cutsite: CutSiteType) -> bool:
       return True

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



if __name__ == "__main__":

    import re

    from Bio.Restriction import RestrictionBatch
    from Bio.Restriction import BamHI, BglII, BsaI, KpnI, Acc65I, XhoII, DpnI, BcgI
    from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

    enzymes = RestrictionBatch((BamHI, BglII, BsaI, KpnI, Acc65I, XhoII))

    # regex_all = re.compile("|".join(enzyme.compsite.pattern for enzyme in enzymes))

    target = "aaGGATCCctAGATCTnn"
    target = "aaAGATCTnnGGATCCwwGGTCTCNNNNNNNN"
    target = "wwGGTCTC12345678aaAGATCTnnGGATCC"
    target = "wwGGTACCGGTCTC12345678aaAGATCTnnGGATCC"
    target = "GGTCTC"
    target = "GAGACC"
    target = "GGATCCnnnGGATCC"

    t = target

    oldobj = oldDseq(t, circular=True)
    oldcs = oldobj.get_cutsites(BamHI)
    oldcsp = oldobj.get_cutsite_pairs(oldcs)
    oldobj.apply_cut(*oldcsp[1])
    oldresult = oldobj.cut(BamHI)

    obj = Dseq(t, circular=True)
    cs = obj.get_cutsites(BamHI)
    csp = obj.get_cutsite_pairs(cs)
    obj.apply_cut(*csp[0])
    result = obj.cut(BamHI)

    # for i in range(len(target) + 1):
    #     t = target[i:]+target[:i]
    #     print(t)
    #     obj = Dseq(t, circular=True)
    #     cs = obj.get_cutsites(BamHI)
    #     csp = obj.get_cutsite_pairs(cs)
    #     result = obj.cut(BamHI)
    #     print()
    #     break

    # for i in range(7):
    #     t = target[i:]+target[:i]
    #     print(t)
    #     print(Dseq(t, circular=True).get_cutsites(BsaI))
    #     print(Dseq(t, circular=True).get_cutsites2(BsaI))
    #     print()

    # m = regex_all.search(target)

    # enz, matched = next((k,v) for k, v in m.groupdict().items() if v)

    # first_enzyme = locals()[enz]

    # enzymes.remove(first_enzyme)

    # regex_wo_first = re.compile("|".join(enzyme.compsite.pattern for enzyme in enzymes))

    # if regex_wo_first.match(target):
    #     raise ValueError("Overlapping enzymes")

    # result = (m.start() + first_enzyme.fst5, m.start() + first_enzyme.fst5 - first_enzyme.ovhg), first_enzyme

    # print(result)


    # Dseq.from_representation("""
    #                               AAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAA
    #                             TTTTTTTTTTTTGCTTTTTTTACGTTTTTTTTTT""").get_cutsites(BcgI)
    # Dseq.from_representation("""
    #                             AAAAAAAAAAAACGAAAAAAATGCAAAAAAAAAAAA
    #                             TTTTTTTTTTTTGCTTTTTTTACGTTTTTTTTTTTT""").get_cutsites(BcgI)
