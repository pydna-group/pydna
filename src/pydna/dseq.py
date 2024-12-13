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
import inspect as _inspect

from pydna.seq import Seq as _Seq
from Bio.Seq import _translate_str, _SeqAbstractBaseClass

from pydna._pretty import pretty_str as _pretty_str
from seguid import ldseguid as _ldseguid
from seguid import cdseguid as _cdseguid

from pydna.utils import rc as _rc
from pydna.utils import flatten as _flatten
from pydna.utils import cuts_overlap as _cuts_overlap
from pydna.utils import complement as _complement

# from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from Bio.Restriction import RestrictionBatch as _RestrictionBatch
from Bio.Restriction import CommOnly
# from itertools import zip_longest

from typing import (
    TYPE_CHECKING,
    List as _List,
    Tuple as _Tuple,
    Union as _Union,
    TypeVar as _TypeVar,
    Iterable as _Iterable,
)

to_watson_table = bytes.maketrans(b"PEXIQFZJpexiqfzj12", b"GATC    gatc     .")
to_crick_table = bytes.maketrans(b"PEXIQFZJpexiqfzj12", b"    GATC    gatc. ")
to_5tail_table = bytes.maketrans(b"GATCgatc", b"QFZJqfzj")
to_3tail_table = bytes.maketrans(b"GATCgatc", b"PEXIpexi")
left_fill_in__table = str.maketrans("PEXIpexi", "GATCgatc")
right_fill_in__table = str.maketrans("QFZJqfzj", "GATCgatc")

#                                       Watson Crick    >   dsIUPAC
ss_to_ds_dct = {(b"P", b"Q"): b"G",   # P      Q            G
                (b"E", b"F"): b"A",   # E      F            A
                (b"X", b"Z"): b"T",   # X      Z            T
                (b"I", b"J"): b"C",   # I      J            C

                (b"p", b"q"): b"g",   # p      q            g
                (b"e", b"f"): b"a",   # e      f            a
                (b"x", b"z"): b"t",   # x      z            t
                (b"i", b"j"): b"c",   # i      j            c

                (b"p", b"Q"): b"g",   # p      Q            g
                (b"e", b"F"): b"a",   # e      F            a
                (b"x", b"Z"): b"t",   # x      Z            t
                (b"i", b"J"): b"c",   # i      J            c

                (b"P", b"q"): b"G",   # P      q            G
                (b"E", b"f"): b"A",   # E      f            A
                (b"X", b"z"): b"T",   # X      z            T
                (b"I", b"j"): b"C",   # I      j            C

                (b"P", b" "): b"P",   # P      q            G
                (b"E", b" "): b"E",   # E      f            A
                (b"X", b" "): b"X",   # X      z            T
                (b"I", b" "): b"I",   # I      j            C

                (b"Q", b"P"): b"G",   # Q      P            G
                (b"F", b"E"): b"A",   # F      E            A
                (b"Z", b"X"): b"T",   # Z      X            T
                (b"J", b"I"): b"C",   # J      I            C

                (b"q", b"p"): b"g",   # q      p            g
                (b"f", b"e"): b"a",   # f      e            a
                (b"z", b"x"): b"t",   # z      x            t
                (b"j", b"i"): b"c",   # j      i            c

                (b"q", b"P"): b"G",   # Q      P            G
                (b"f", b"E"): b"A",   # F      E            A
                (b"z", b"X"): b"T",   # Z      X            T
                (b"j", b"I"): b"C",   # J      I            C

                (b"Q", b"p"): b"g",   # q      p            g
                (b"F", b"e"): b"a",   # f      e            a
                (b"Z", b"x"): b"t",   # z      x            t
                (b"J", b"i"): b"c",   # j      i            c

                (b" ", b"Q"): b"Q",   # Q      P            G
                (b" ", b"F"): b"F",   # F      E            A
                (b" ", b"Z"): b"Z",   # Z      X            T
                (b" ", b"J"): b"J",   # J      I            C

                (b"G", b" "): b"P",
                (b"A", b" "): b"E",
                (b"T", b" "): b"X",
                (b"C", b" "): b"I",
                (b"g", b" "): b"p",
                (b"a", b" "): b"e",
                (b"t", b" "): b"x",
                (b"c", b" "): b"i",

                (b" ", b"G"): b"J",
                (b" ", b"A"): b"Z",
                (b" ", b"T"): b"F",
                (b" ", b"C"): b"Q",
                (b" ", b"g"): b"j",
                (b" ", b"a"): b"z",
                (b" ", b"t"): b"f",
                (b" ", b"c"): b"q",

                (b"G", b"C"): b"G",
                (b"A", b"T"): b"A",
                (b"T", b"A"): b"T",
                (b"C", b"G"): b"C",
                (b"g", b"c"): b"g",
                (b"a", b"t"): b"a",
                (b"t", b"a"): b"t",
                (b"c", b"g"): b"c", }

if TYPE_CHECKING:
    from Bio.Restriction import AbstractCut as _AbstractCut


# To represent any subclass of Dseq
DseqType = _TypeVar("DseqType", bound="Dseq")
EnzymesType = _TypeVar("EnzymesType", _RestrictionBatch, _Iterable["_AbstractCut"], "_AbstractCut")
CutSiteType = _Tuple[_Tuple[int, int], _Union["_AbstractCut", None]]


class Dseq(_Seq):
    """docstring.

    here
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
        self.circular = circular
        self.pos = pos

    @classmethod
    def from_representation(cls, dsdna: str, *args, **kwargs):
        obj = cls.__new__(cls)  # Does not call __init__
        obj.circular = False
        obj.pos = 0
        dsdna = _inspect.cleandoc(dsdna)
        # breakpoint()
        watson, crick = [ln for ln in dsdna.splitlines() if ln.strip() and not ln.strip().startswith("Dseq(")]
        watson = f"{watson:<{len(crick)}}"
        crick = f"{crick:<{len(watson)}}"
        data = bytearray()
        # breakpoint()
        for w, c in zip(watson, crick):
            # print(w,c, ss_to_ds_dct[w.encode("ascii"), c.encode("ascii")])
            data.extend(ss_to_ds_dct[w.encode("ascii"), c.encode("ascii")])
        obj._data = data
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
            data.extend(ss_to_ds_dct[w.encode("ascii"), c.encode("ascii")])

        return cls(data)

    @property
    def watson(self):
        return self._data.translate(to_watson_table).strip().decode("ascii")

    @property
    def crick(self):
        return _complement(self._data.translate(to_crick_table).strip().decode("ascii"))[::-1]

    @property
    def ovhg(self):
        ohw = len(_re.match(b"^[PEXIpexi]*",self._data).group(0))
        ohc = len(_re.match(b"^[QFZJqfzj]*",self._data).group(0))
        return -ohw or ohc

    def __getitem__(self, sl: slice) -> "Dseq":
        """Returns a subsequence. This method is used by the slice notation"""
        return self.__class__(super().__getitem__(sl))

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

    def __eq__(self, other: DseqType) -> bool:
        """Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties. This comparison is case
        insensitive.

        """
        try:
            same = self._data.lower() == other._data.lower() and self.circular == other.circular
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
        m = _re.match(b"([PEXIpexi]*)([QFZJqfzj]*)(?=[GATCgatc])(.*)(?<=[GATCgatc])([PEXIpexi]*)([QFZJqfzj]*)|([PEXIpexiQFZJqfzj]+)", self._data)
        result = m.groups() if m else (b"",) * 7
        sticky_left5, sticky_left3, middle, sticky_right5, sticky_right3, single = result
        if len(self) > self.trunc:
            sticky_left5 = sticky_left5[:4] + b"22" + sticky_left5[-4:] if sticky_left5 and len(sticky_left5) > 10 else sticky_left5
            sticky_left3 = sticky_left3[:4] + b"11" + sticky_left3[-4:] if sticky_left3 and len(sticky_left3) > 10 else sticky_left3
            middle = middle[:4] + b".." + middle[-4:] if middle and len(middle) > 30 else middle
            sticky_right5 = sticky_right5[:4] + b"22" + sticky_right5[-4:] if sticky_right5 and len(sticky_right5) > 10 else sticky_right5
            sticky_right3 = sticky_right3[:4] + b"11" + sticky_right3[-4:] if sticky_right3 and len(sticky_right3) > 10 else sticky_right3
        r = (sticky_left5 or sticky_left3 or b"") + (middle or b"") + (sticky_right5 or sticky_right3 or single or b"")
        return _pretty_str(f"{header}\n{r.translate(to_watson_table).decode().rstrip()}\n{_complement(r.translate(to_crick_table)).decode()}")

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

        ms = _re.fullmatch(b"(?:.*)([GATCgatc])([QFZJqfzj]+|[PEXIpexi]+)", self._data)
        mo = _re.fullmatch(b"([PEXIpexi]+|[QFZJqfzj]+)([GATCgatc])(?:.*)", self._data)

        sticky_left = ms.group(2) if ms else b""
        sticky_right = mo.group(1) if mo else b""

        sticky_left_just = bytes(f"{sticky_left.decode():<{len(sticky_right)}}", 'utf8')
        sticky_right_just = bytes(f"{sticky_right.decode():>{len(sticky_left)}}", 'utf8')

        assert len(sticky_left_just) == len(sticky_right_just)

        junction = b"".join([ss_to_ds_dct.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_left_just, sticky_right_just)])

        if b"-" in junction:
            raise TypeError("DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!")

        return self.__class__(junction + self._data[len(sticky_left) or None:-len(sticky_right) or None], circular=True)

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

    def __add__(self: DseqType, other: DseqType) -> DseqType:

        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        ms = _re.fullmatch(b"(?:.*)([GATCgatc])([QFZJqfzj]+|[PEXIpexi]+)", self._data)
        mo = _re.fullmatch(b"([PEXIpexi]+|[QFZJqfzj]+)([GATCgatc])(?:.*)", other._data)

        sticky_self = ms.group(2) if ms else b""
        sticky_other = mo.group(1) if mo else b""

        sticky_self_just = bytes(f"{sticky_self.decode():<{len(sticky_other)}}", 'utf8')
        sticky_other_just = bytes(f"{sticky_other.decode():>{len(sticky_self)}}", 'utf8')

        assert len(sticky_self_just) == len(sticky_other_just)

        junction = b"".join([ss_to_ds_dct.get((bytes([w]), bytes([c])), b"-") for w, c in zip(sticky_self_just, sticky_other_just)])

        if b"-" in junction:
            raise TypeError("sticky ends not compatible!")

        return self.__class__(self._data[:-len(sticky_self) or None] + junction + other._data[len(sticky_other) or None:])

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

        """
        mleft = _re.fullmatch(b"([PEXIpexi]+)([GATCgatc])(?:.*)", self._data)
        sticky_left = mleft.group(1).decode("ASCII")[::-1] if mleft else ""

        allowed = _complement(nucleotides.upper())
        newletters = list(sticky_left)

        for i, letter in enumerate(sticky_left.translate(left_fill_in__table)):
            if letter.upper() in allowed:
                newletters[i] = letter
            else:
                break

        end = ''.join(newletters[::-1])

        return self.__class__(end.encode("ASCII") + self._data[len(sticky_left) or None:], circular=False)

    def _fill_in_five_prime(self: DseqType, nucleotides: str) -> str:
        """
        NNNNQFZJ

        NNNN
        NNNNCTAG

        """

        mright = _re.fullmatch(b"(?:.*)([GATCgatc])([QFZJqfzj]+)", self._data)
        sticky_right = mright.group(2).decode("ASCII") if mright else ""

        allowed = nucleotides.upper()
        newletters = list(sticky_right)

        for i, letter in enumerate(sticky_right.translate(right_fill_in__table)):
            if letter.upper() in allowed:
                newletters[i] = letter
            else:
                break

        end = ''.join(newletters)

        return self.__class__(self._data[:-len(sticky_right) or None] + end.encode("ASCII"), circular=False)

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
            nucleotides = "GATCgatc"

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

        return self.__class__(self._data.strip(b"PEXIQFZJpexiqfzj"), circular=False)

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

        nucleotides = nucleotides if isinstance(nucleotides, bytes) else nucleotides.encode('ascii')

        strp = self.lstrip(b"QFZJqfzj").rstrip(b"PEXIpexi").fill_in(nucleotides.decode("ascii"))  # remove 3' sticky ends on both sides

        to_remove = bytes(set(b"GATCgatc") - set(nucleotides.upper()+nucleotides.lower()))

        def repl(m):
            breakpoint()
            s2 = m.group(2).translate(to_3tail_table) or b""
            s4 = m.group(4).translate(to_5tail_table) or b""
            sub = m.group(1) + s2 + m.group(3) + s4 + m.group(5)
            print("hej")
            return sub

      # new = _re.sub(b"([PEXIpexi]*)(?=[GATCgatc])([%b]+)(.+?)([%b]+)(?<=[GATCgatc])([QFZJqfzj]*)"% (_complement(to_remove), to_remove), repl, strp._data)
        new = _re.sub(b"([PEXIpexi]*)(?=[GATCgatc])([%b]+)(.*?)([%b]*)(?<=[GATCgatc])([QFZJqfzj]*)"% (_complement(to_remove), to_remove), repl, strp._data)
        return self.__class__(new).fill_in(nucleotides.decode("ascii"))

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
            w = self._data.translate(to_watson_table).strip().decode("ascii").upper().replace(" ", "-")
            c = _complement(self._data.translate(to_crick_table).strip().decode("ascii")).upper()[::-1].replace(" ", "-")
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
        return (self._data[0] not in b"PEXIQFZJpexiqfzj") and (self._data[-1] not in b"PEXIQFZJpexiqfzj") and not self.circular

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
        cutsite_pairs = [((0, 0), None)] + self.get_cutsites(*enzymes) + [((len(self), len(self)), None)]

        frags = []

        for ((w1, c1), e), ((w2, c2), e) in _itertools.pairwise(cutsite_pairs):

            table1 = {1.0: to_5tail_table, -1.0: to_3tail_table}[_math.copysign(1, w1 - c1)]
            table2 = {1.0: to_3tail_table, -1.0: to_5tail_table}[_math.copysign(1, w2 - c2)]

            frags.append(self._data[min((w1, c1)):max((w1, c1))].translate(table1)
                         + self._data[max((w1, c1)):min((w2, c2))]
                         + self._data[min((w2, c2)):max((w2, c2))].translate(table2))

        return tuple(self.__class__(frag) for frag in frags)

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
            cuts = [((i.start() + e.fst5, i.start() + e.size + e.fst3), e) for i in _re.finditer(e.compsite, str(self))]
            out.extend(cuts)

        return sorted([cutsite for cutsite in out])

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
    # import os as _os

    # cached = _os.getenv("pydna_cached_funcs", "")
    # _os.environ["pydna_cached_funcs"] = ""
    # import doctest

    # doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    # _os.environ["pydna_cached_funcs"] = cached

    # for x, y in [("abc", "abc"),("abc", "ab"),("ab", "abc")]:

    #     a = f"{x:<{len(y)}}"
    #     b = f"{y:<{len(x)}}"
    #     print(repr(a))
    #     print(repr(b))
    #     print()

    fiveoh = Dseq('PEXIaaaQFZJ')
    assert str(fiveoh + fiveoh) == 'PEXIaaaGATCaaaQFZJ'
    threeoh = Dseq('QFZJQtttPEXIP')
    assert str(threeoh + threeoh) == "QFZJQtttGATCGtttPEXIP"

    assert repr(Dseq("AIXEP") + Dseq("JZFQA")) == "Dseq(-6)\nACTAGA\nTGATCT"
    assert repr(Dseq("AIXEP") + Dseq("ZFQA")) == "Dseq(-6)\nACTAGA\nT ATCT"
    assert repr(Dseq("AIXE") + Dseq("JZFQA")) == "Dseq(-6)\nACTA A\nTGATCT"
    assert repr(Dseq("AP") + Dseq("QA")) == "Dseq(-3)\nAGA\nTCT"
    assert repr(Dseq("AE") + Dseq("FA")) == "Dseq(-3)\nAAA\nTTT"
    assert repr(Dseq("AX") + Dseq("ZA")) == "Dseq(-3)\nATA\nTAT"
    assert repr(Dseq("AI") + Dseq("JA")) == "Dseq(-3)\nACA\nTGT"
    assert repr(Dseq("APP") + Dseq("QA")) == "Dseq(-4)\nAGGA\nT CT"
    assert repr(Dseq("AEE") + Dseq("FA")) == "Dseq(-4)\nAAAA\nT TT"
    assert repr(Dseq("AXX") + Dseq("ZA")) == "Dseq(-4)\nATTA\nT AT"
    assert repr(Dseq("AII") + Dseq("JA")) == "Dseq(-4)\nACCA\nT GT"

    assert repr(Dseq("AP") + Dseq("QQA")) == "Dseq(-4)\nAG A\nTCCT"
    assert repr(Dseq("AE") + Dseq("FFA")) == "Dseq(-4)\nAA A\nTTTT"
    assert repr(Dseq("AX") + Dseq("ZZA")) == "Dseq(-4)\nAT A\nTAAT"
    assert repr(Dseq("AI") + Dseq("JJA")) == "Dseq(-4)\nAC A\nTGGT"

    from Bio.Restriction import XmaI, SmaI, KpnI, Acc65I

    s = Dseq("CCCGGGGCATCGTAGTGATCGGTACC") #  blunt

    assert s.looped()._data == s._data

    a,b,c = s.cut([XmaI, Acc65I])

    assert repr(a+b+c) == repr(s)

    a,b,c = s.cut([SmaI, Acc65I])

    assert repr(a+b+c) == repr(s)

    a,b,c = s.cut([SmaI, KpnI])

    assert repr(a+b+c) == repr(s)

    s = Dseq("PPCCCGGGGCATCGTAGTGATCGGTACC")

    a,b,c = s.cut([XmaI, Acc65I])

    assert repr(a+b+c) == repr(s)

    s = Dseq("QQCCCGGGGCATCGTAGTGATCGGTACC")

    a,b,c = s.cut([XmaI, Acc65I])

    assert repr(a+b+c) == repr(s)

    s = Dseq("PEXIAAAQFZJ")

    assert s.looped()._data == b'GATCAAA'

    s = Dseq("PEXIAAAQFZJ")

    assert s._fill_in_three_prime("gatc")._data == b'GATCAAAQFZJ'
    assert s._fill_in_three_prime("gat")._data == b'PATCAAAQFZJ'
    assert s._fill_in_three_prime("ga")._data == b'PETCAAAQFZJ'
    assert s._fill_in_three_prime("g")._data == b'PEXCAAAQFZJ'
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_three_prime("atc")._data == s._data
    assert s._fill_in_three_prime("at")._data == s._data
    assert s._fill_in_three_prime("a")._data == s._data
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_five_prime("gatc")._data == b'PEXIAAAGATC'
    assert s._fill_in_five_prime("gat")._data ==  b'PEXIAAAGATJ'
    assert s._fill_in_five_prime("ga")._data ==   b'PEXIAAAGAZJ'
    assert s._fill_in_five_prime("g")._data ==    b'PEXIAAAGFZJ'
    assert s._fill_in_five_prime("")._data == s._data

    assert s._fill_in_five_prime("atc")._data == s._data
    assert s._fill_in_five_prime("at")._data == s._data
    assert s._fill_in_five_prime("a")._data == s._data
    assert s._fill_in_five_prime("")._data == s._data

    assert s.fill_in("gatc")._data == b'GATCAAAGATC'
    assert s.fill_in("gat")._data == b'PATCAAAGATJ'
    assert s.fill_in("ga")._data == b'PETCAAAGAZJ'
    assert s.fill_in("g")._data == b'PEXCAAAGFZJ'
    assert s.fill_in("")._data == s._data

    assert s == Dseq.from_representation(repr(s))

    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=2)._data == Dseq("FFAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=2)._data == Dseq("EEAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=-2)._data == Dseq("FFAAFF")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=-2)._data == Dseq("EEAAFF")._data

    assert repr(Dseq("XIpexiPEXIpexiAqfzjQFZJqfzjQF")) == ("Dseq(-29)\n"
                                                           "TCgatcGATCgatcA\n"
                                                           "              TctagCTAGctagCT")

    assert repr(Dseq("XIpexiPEXIpexiAAqfzjQFZJqfzjQF")) == ("Dseq(-30)\n"
                                                             "TCgatcGATCgatcAA\n"
                                                             "              TTctagCTAGctagCT")

    assert repr(Dseq("EXIpexiPEXIpexiAqfzjQFZJqfzjQFZ")) == ("Dseq(-31)\n"
                                                             "ATCg..gatcA\n"
                                                             "          Tctag..gCTA")

    assert repr(Dseq("XIpexiPEXIpexiAaAqfzjQFZJqfzjQF")) == ("Dseq(-31)\n"
                                                             "TCga..gatcAaA\n"
                                                             "          TtTctag..agCT")



    dsdna = """
               Dseq(
                GGATCC
               aCCTAGGg"""


    assert Dseq.from_representation(dsdna)._data == b"zGGATCCj"

    dsdna = """Dseq(
               aGGATCCg
                CCTAGG"""

    assert Dseq.from_representation(dsdna)._data == b"eGGATCCp"

    self = Dseq('PEXIGGATCCQFZJ')
    # self = Dseq('PAGAJ')
    # self = Dseq('QFZJGGATCCPEXI')


    Dseq("pexiAqfzj").T4(b"GATC")
