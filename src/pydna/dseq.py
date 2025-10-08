#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Provides the Dseq class for handling double stranded DNA sequences.

Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class
which can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology.
"""


import copy as _copy
import re as _re
import sys as _sys
import math as _math
import inspect as _inspect
from typing import List as _List, Tuple as _Tuple, Union as _Union

from Bio.Restriction import RestrictionBatch as _RestrictionBatch
from Bio.Restriction import CommOnly

from seguid import ldseguid as _ldseguid
from seguid import cdseguid as _cdseguid

from pydna.seq import Seq as _Seq
from Bio.Seq import _translate_str, _SeqAbstractBaseClass
from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import rc as _rc
from pydna.utils import flatten as _flatten
from pydna.utils import cuts_overlap as _cuts_overlap

from pydna.alphabet import basepair_dict
from pydna.alphabet import dscode_to_watson_table
from pydna.alphabet import dscode_to_crick_table

# from pydna.alphabet import watson_tail_to_dscode_table
# from pydna.alphabet import crick_tail_to_dscode_table

from pydna.alphabet import dscode_to_full_sequence_table
from pydna.alphabet import dscode_to_watson_tail_table
from pydna.alphabet import dscode_to_crick_tail_table

from pydna.alphabet import get_parts
from pydna.alphabet import representation_tuple
from pydna.alphabet import dsbreaks

from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from pydna.types import DseqType, EnzymesType, CutSiteType


# Sequences larger than this gets a truncated representation.
length_limit_for_repr = 30


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


class Dseq(_Seq):
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
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/dist-packages/pydna_/dsdna.py", line 169, in __init__
        else:
    ValueError: ovhg defined without crick strand!

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
    >>> Dseq("aaa", "ttt", ovhg=1, circular = True)
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
        watson: _Union[str, bytes],
        crick: _Union[str, bytes, None] = None,
        ovhg=None,
        circular=False,
        pos=0,
    ):
        if isinstance(watson, bytes):
            # watson is decoded to a string if needed.
            watson = watson.decode("ascii")
        if isinstance(crick, bytes):
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
                limit = int(_math.log(len(watson)) / _math.log(4))
                olaps = _common_sub_strings(
                    str(watson).lower(),
                    str(_rc(crick).lower()),
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
                This means that annealing is ambigous. USer should provide
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
                    raise ValueError()
            data = "".join(data)
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
    def from_string(
        cls,
        dna: str,
        *args,
        circular=False,
        **kwargs,
    ):
        obj = cls.__new__(cls)
        obj.circular = circular
        obj.pos = 0
        obj._data = dna.encode("ascii")
        return obj

    @classmethod
    def from_representation(cls, dsdna: str, *args, **kwargs):
        obj = cls.__new__(cls)
        obj.circular = False
        obj.pos = 0
        clean = _inspect.cleandoc("\n" + dsdna)
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

        # watson = f"{watson:<{len(crick)}}"
        # crick = f"{crick:<{len(watson)}}"
        # data = bytearray()
        # for w, c in zip(watson, crick):
        #     try:
        #         data.extend(basepair_dict[w, c])
        #     except KeyError as err:
        #         print(f"Base mismatch in representation {err}")
        #         raise
        # obj._data = bytes(data)
        # return obj

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
    def ovhg(self) -> int:
        """
        The 5' overhang of the lower strand compared the the upper.

        See module docstring for more information.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        parts = get_parts(self._data.decode("ascii"))
        if parts.single_watson or parts.single_crick:
            return None
        return -len(parts.sticky_left5) or len(parts.sticky_left3)

    def to_blunt_string(self):
        """
        The sequence as a blunt ended string.


        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self._data.decode("ascii").translate(dscode_to_full_sequence_table)

    __str__ = to_blunt_string

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

    def find(
        self, sub: _Union[_SeqAbstractBaseClass, str, bytes], start=0, end=_sys.maxsize
    ) -> int:
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

    def __contains__(self, item: [str, bytes]):
        if isinstance(item, bytes):
            item = item.decode("ascii")
        return item in self.to_blunt_string()

    def __getitem__(self, sl: [slice, int]) -> "Dseq":
        """Method is used by the slice notation"""
        if isinstance(sl, int):
            # slice(start, stop, step) where stop = start+1
            sl = slice(sl, sl + 1, 1)
        sl = slice(sl.start or 0, sl.stop or len(self), sl.step)
        if self.circular:
            if sl.start is None and sl.stop is None:
                # Returns a linear copy using slice obj[:]
                return self.quick(self._data[sl])
            elif sl.start == sl.stop:
                # Returns a shifted object
                shift = sl.start % len(self)
                return self.quick(self._data[shift:] + self._data[:shift])
            elif (
                sl.start > sl.stop
                and 0 <= sl.start <= len(self)
                and 0 <= sl.stop <= len(self)
            ):
                # Returns the circular slice spanning the origin
                return self.quick(self._data[sl.start :] + self._data[: sl.stop])
        return super().__getitem__(sl)

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

    def __repr__(self, lim: int = length_limit_for_repr):

        header = f"{self.__class__.__name__}({({False: '-', True: 'o'}[self.circular])}{len(self)})"

        w, c = representation_tuple(
            self._data.decode("ascii"), length_limit_for_repr=length_limit_for_repr
        )

        return _pretty_str(header + "\n" + w + "\n" + c)

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
        return Dseq.quick(_rc(self._data), circular=self.circular)

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

        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()
        if type5 == type3 and str(sticky5) == str(_rc(sticky3)):
            return self.__class__(self.watson, circular=True)
        else:
            raise TypeError(
                "DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!"
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
            "tolinear method is obsolete; "
            "please use obj[:] "
            "instead of obj.tolinear().",
            _PydnaDeprecationWarning,
        )
        if not self.circular:
            raise TypeError("DNA is not circular.\n")
        selfcopy = _copy.deepcopy(self)
        selfcopy.circular = False
        return selfcopy

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
        parts = get_parts(self._data.decode("ascii"))

        sticky5 = parts.sticky_left5.translate(dscode_to_watson_table)

        sticky3 = parts.sticky_left3.translate(dscode_to_crick_table)[::-1]

        single_watson = parts.single_watson.translate(dscode_to_watson_table)

        single_crick = parts.single_crick.translate(dscode_to_crick_table)[::-1]

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

        parts = get_parts(self._data.decode("ascii"))

        sticky5 = parts.sticky_right5.translate(dscode_to_crick_table)[::-1]

        sticky3 = parts.sticky_right3.translate(dscode_to_watson_table)

        single_watson = parts.single_watson.translate(dscode_to_watson_table)

        single_crick = parts.single_crick.translate(dscode_to_crick_table)[::-1]

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

    def watson_ovhg(self) -> int:
        """Overhang of the watson strand at the right side."""
        parts = get_parts(self._data.decode("ascii"))
        if parts.single_watson or parts.single_crick:
            return None
        return -len(parts.sticky_right5) or len(parts.sticky_right3)

    def __add__(self: DseqType, other: DseqType) -> DseqType:

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

        self_type, self_tail = self.three_prime_end()
        other_type, other_tail = other.five_prime_end()

        if self_type == other_type and str(self_tail) == str(_rc(other_tail)):
            return self.__class__(
                self.watson + other.watson, other.crick + self.crick, self.ovhg
            )
        raise TypeError("sticky ends not compatible!")

    def __mul__(self: DseqType, number: int) -> DseqType:
        if not isinstance(number, int):
            raise TypeError(
                "TypeError: can't multiply Dseq by non-int of type {}".format(
                    type(number)
                )
            )
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

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ) -> _Seq:
        return _Seq(
            _translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap)
        )

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
        return Dseq(
            self.watson[
                max(0, -self.ovhg) : min(len(self.watson), len(self.crick) - self.ovhg)
            ]
        )

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

    def no_cutters(
        self, batch: _Union[_RestrictionBatch, None] = None
    ) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch not cutting sequence."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if not sitelist}
        return _RestrictionBatch(ncut)

    def unique_cutters(
        self, batch: _Union[_RestrictionBatch, None] = None
    ) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence once."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=1, batch=batch)

    once_cutters = unique_cutters  # alias for unique_cutters

    def twice_cutters(
        self, batch: _Union[_RestrictionBatch, None] = None
    ) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence twice."""
        if batch is None:
            batch = CommOnly
        return self.n_cutters(n=2, batch=batch)

    def n_cutters(
        self, n=3, batch: _Union[_RestrictionBatch, None] = None
    ) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting n times."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if len(sitelist) == n}
        return _RestrictionBatch(ncut)

    def cutters(
        self, batch: _Union[_RestrictionBatch, None] = None
    ) -> _RestrictionBatch:
        """Enzymes in a RestrictionBatch cutting sequence at least once."""
        if batch is None:
            batch = CommOnly
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if sitelist}
        return _RestrictionBatch(ncut)

    def seguid(self) -> str:
        """SEGUID checksum for the sequence."""
        if self.circular:
            cs = _cdseguid(
                self.watson.upper(), self.crick.upper(), alphabet="{DNA-extended}"
            )
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
        parts = get_parts(self._data.decode("ascii"))

        return not any(
            (
                parts.sticky_right5,
                parts.sticky_right3,
                parts.sticky_left3,
                parts.sticky_left5,
                self.circular,
            )
        )

    def cas9(self, RNA: str) -> _Tuple[slice, ...]:
        """docstring."""
        bRNA = RNA.encode("ascii")
        slices = []
        cuts = [0]
        for m in _re.finditer(bRNA, self._data):
            cuts.append(m.start() + 17)
        cuts.append(len(self))
        slices = tuple(slice(x, y, 1) for x, y in zip(cuts, cuts[1:]))
        return slices

    def terminal_transferase(self, nucleotides="a") -> "Dseq":
        """docstring."""
        ovhg = self.ovhg
        if self.ovhg >= 0:
            ovhg += len(nucleotides)
        return Dseq(self.watson + nucleotides, self.crick + nucleotides, ovhg)

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
        if (
            len(recognition_site) == 0
            or recognition_site.ovhg != 0
            or recognition_site.watson_ovhg() != 0
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
                    or recognition_site.watson_ovhg() != 0
                ):
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

    def get_ss_meltsites(self: DseqType, length) -> _List[CutSiteType]:

        if length < 1:
            return ()

        regex = (
            f"(?P<watson>((?<=[PEXIpexi]))([GATCgatcUuOo]{{1,{length}}})((?=[^QFZJqfzjGATCgatcUuOo])))|"
            f"(?P<crick>((?<=[QFZJqfzj]))([GATCgatcUuOo]{{1,{length}}})((?=[^PEXIpexiGATCgatcUuOo])))"
        )

        regex = _re.compile(regex.encode("ascii"))

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

        regex = (
            f"(?P<watson>((?<=[PEXIpexi])|^)([GATCgatcUuOo]{{1,{length}}})((?=[^PEXIpexiGATCgatcUuOo])|$))|"
            f"(?P<crick>((?<=[QFZJqfzj])|^)([GATCgatcUuOo]{{1,{length}}})((?=[^QFZJqfzjGATCgatcUuOo])|$))"
        )

        regex = _re.compile(regex.encode("ascii"))

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
        NNNN----
        NNNNCTAG

        NNNNGATC
        NNNNCTAG

        NNNNGATC
        NNNN----

        NNNNGATC
        NNNNCTAG

        """

        p = get_parts(self._data.decode("ascii"))

        ds_stuffer = (p.sticky_right5 or p.sticky_right3).translate(
            dscode_to_full_sequence_table
        )

        result = (p.sticky_left5 or p.sticky_left3) + p.middle + ds_stuffer

        return self.__class__(result, circular=False)

    def cast_to_ds_left(self):
        """
        GATCNNNN
            NNNN

        GATCNNNN
        CTAGNNNN

            NNNN
        CTAGNNNN

        GATCNNNN
        CTAGNNNN

        """

        p = get_parts(self._data.decode("ascii"))

        ds_stuffer = (p.sticky_left5 or p.sticky_left3).translate(
            dscode_to_full_sequence_table
        )

        result = ds_stuffer + p.middle + (p.sticky_right5 or p.sticky_right3)

        return self.__class__(result, circular=False)

    def get_cut_parameters(
        self, cut: _Union[CutSiteType, None], is_left: bool
    ) -> _Tuple[int, int, int]:
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

    def melt(self, length):
        if not length or length < 1:
            return tuple()

        new, strands = self.shed_ss_dna(length)

        cutsites = new.get_ds_meltsites(length)

        cutsite_pairs = self.get_cutsite_pairs(cutsites)

        result = tuple(new.apply_cut(*cutsite_pair) for cutsite_pair in cutsite_pairs)

        result = tuple([new]) if strands and not result else result

        return strands + result

    def shed_ss_dna(self, length):

        new, strands, intervals = self._shed_ss_dna(length)

        return new, strands

    def _shed_ss_dna(self, length):

        watsonnicks, cricknicks = self.get_ss_meltsites(length)

        watsonstrands = []
        crickstrands = []

        new = self._data

        for x, y in watsonnicks:
            stuffer = new[x:y]
            ss = Dseq.quick(stuffer.translate(dscode_to_watson_tail_table))
            new = new[:x] + stuffer.translate(dscode_to_crick_tail_table) + new[y:]
            watsonstrands.append((x, y, ss))

        for x, y in cricknicks:
            stuffer = new[x:y]
            ss = Dseq.quick(stuffer.translate(dscode_to_crick_tail_table))
            new = new[:x] + stuffer.translate(dscode_to_watson_tail_table) + new[y:]
            crickstrands.append((x, y, ss))

        ordered_strands = sorted(watsonstrands + crickstrands)

        strands = []

        for x, y, ss in ordered_strands:
            seq = (
                ss._data[::-1].translate(dscode_to_watson_table).strip()
                or ss._data.translate(dscode_to_crick_table).strip()
            )
            strands.append(_Seq(seq))

        intervals = tuple((x, y) for x, y, ss in ordered_strands)

        return Dseq(new), tuple(strands), intervals

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
            str(
                self.reverse_complement()[
                    len(self) - right_crick : len(self) - left_crick
                ]
            ),
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
