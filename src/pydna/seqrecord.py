#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""
A subclass of the Biopython SeqRecord class.

Has a number of extra methods and uses
the :class:`pydna._pretty_str.pretty_str` class instread of str for a
nicer output in the IPython shell.
"""


from Bio.SeqFeature import SeqFeature as _SeqFeature
from pydna._pretty import pretty_str as _pretty_str

from pydna.seq import ProteinSeq as _ProteinSeq
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings

from Bio.Data.CodonTable import TranslationError as _TranslationError
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from pydna.seq import Seq as _Seq
from pydna._pretty import PrettyTable as _PrettyTable

import re as _re
import pickle as _pickle
from copy import copy as _copy

from pydna import _PydnaWarning
from warnings import warn as _warn

import logging as _logging
import datetime

_module_logger = _logging.getLogger("pydna." + __name__)


class SeqRecord(_SeqRecord):
    """
    A subclass of the Biopython SeqRecord class.

    Has a number of extra methods and uses
    the :class:`pydna._pretty_str.pretty_str` class instread of str for a
    nicer output in the IPython shell.
    """

    def __init__(
        self,
        seq,
        id="id",
        name="name",
        description="description",
        dbxrefs=None,
        features=None,
        annotations=None,
        letter_annotations=None,
    ):
        if isinstance(seq, str):
            seq = _Seq(seq)
        super().__init__(
            seq,
            id=id,
            name=name,
            description=description,
            dbxrefs=dbxrefs,
            features=features,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )
        self._fix_attributes()

    def _fix_attributes(self):
        self.id = _pretty_str(self.id)
        self.name = _pretty_str(self.name)
        self.description = _pretty_str(self.description)

        self.annotations.update({"molecule_type": "DNA"})
        self.map_target = None

        if not hasattr(self.seq, "transcribe"):
            self.seq = _Seq(self.seq)

        self.seq._data = b"".join(self.seq._data.split())  # remove whitespaces
        self.annotations = {_pretty_str(k): _pretty_str(v) for k, v in self.annotations.items()}

    @classmethod
    def from_Bio_SeqRecord(clc, sr: _SeqRecord):
        """Creates a pydnaSeqRecord from a Biopython SeqRecord."""
        # https://stackoverflow.com/questions/15404256/changing-the-\
        # class-of-a-python-object-casting
        sr.__class__ = clc
        sr._fix_attributes()
        return sr

    @property
    def locus(self):
        """Alias for name property."""
        return self.name

    @locus.setter
    def locus(self, value):
        """Alias for name property."""
        if len(value) > 16:
            shortvalue = value[:16]
            _warn(
                ("locus property {} truncated" "to 16 chars {}").format(value, shortvalue),
                _PydnaWarning,
                stacklevel=2,
            )
            value = shortvalue
        self.name = value
        return

    @property
    def accession(self):
        """Alias for id property."""
        return self.id

    @accession.setter
    def accession(self, value):
        """Alias for id property."""
        self.id = value
        return

    @property
    def definition(self):
        """Alias for description property."""
        return self.description

    @definition.setter
    def definition(self, value):
        """Alias for id property."""
        self.description = value
        return

    def reverse_complement(self, *args, **kwargs):
        """Return the reverse complement of the sequence."""
        answer = super().reverse_complement(*args, **kwargs)
        answer = type(self).from_Bio_SeqRecord(answer)
        return answer

    rc = reverse_complement

    def isorf(self, table=1):
        """Detect if sequence is an open reading frame (orf) in the 5'-3'.

        direction.

        Translation tables are numbers according to the NCBI numbering [#]_.

        Parameters
        ----------
        table  : int
            Sets the translation table, default is 1 (standard code)

        Returns
        -------
        bool
            True if sequence is an orf, False otherwise.


        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.isorf()
        True
        >>> b=SeqRecord("atgaaa")
        >>> b.isorf()
        False
        >>> c=SeqRecord("atttaa")
        >>> c.isorf()
        False

        References
        ----------
        .. [#] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

        """
        try:
            self.seq.translate(table=table, cds=True)
        except _TranslationError:
            return False
        else:
            return True

    def translate(self):
        """docstring."""
        p = super().translate()
        return ProteinSeqRecord(_ProteinSeq(p.seq[:-1]))

    def add_colors_to_features_for_ape(self):
        """Assign colors to features.

        compatible with
        the `ApE editor <http://jorgensen.biology.utah.edu/wayned/ape/>`_.

        """
        cols = (
            "#66ffa3",
            "#84ff66",
            "#e0ff66",
            "#ffc166",
            "#ff6666",
            "#ff99d6",
            "#ea99ff",
            "#ad99ff",
            "#99c1ff",
            "#99ffff",
            "#99ffc1",
            "#adff99",
            "#eaff99",
            "#ffd699",
            "#ff9999",
            "#ffccea",
            "#f4ccff",
            "#d6ccff",
            "#cce0ff",
            "#ccffff",
            "#ccffe0",
            "#d6ffcc",
            "#f4ffcc",
            "#ffeacc",
            "#ffcccc",
            "#ff66c1",
            "#e066ff",
            "#8466ff",
            "#66a3ff",
            "#66ffff",
        )

        for i, f in enumerate(self.features):
            f.qualifiers["ApEinfo_fwdcolor"] = [cols[i % len(cols)]]
            f.qualifiers["ApEinfo_revcolor"] = [cols[::-1][i % len(cols)]]

    def add_feature(self, x=None, y=None, seq=None, type_="misc", strand=1, *args, **kwargs):
        """Add a feature of type misc to the feature list of the sequence.

        Parameters
        ----------
        x  : int
            Indicates start of the feature
        y  : int
            Indicates end of the feature

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.features
        []
        >>> a.add_feature(2,4)
        >>> a.features
        [SeqFeature(SimpleLocation(ExactPosition(2),
                                   ExactPosition(4),
                                   strand=1),
                    type='misc',
                    qualifiers=...)]
        """
        qualifiers = {}
        qualifiers.update(kwargs)

        if seq:
            if hasattr(seq, "seq"):
                seq = seq.seq
                if hasattr(seq, "watson"):
                    seq = str(seq.watson).lower()
                else:
                    seq = str(seq).lower()
            else:
                seq = str(seq).lower()
            x = self.seq.lower().find(seq)
            if x == -1:
                raise TypeError("Could not find {}".format(seq))
            y = x + len(seq)
        else:
            x = x or 0
            y = y or len(self)

        if "label" not in qualifiers:
            qualifiers["label"] = ["ft{}".format(y - x)]

            if self[x:y].isorf() or self[x:y].reverse_complement().isorf():
                qualifiers["label"] = ["orf{}".format(y - x)]

        try:
            location = _SimpleLocation(x, y, strand=strand)
        except ValueError as err:
            if self.circular:
                location = _CompoundLocation(
                    (
                        _SimpleLocation(x, self.seq.length, strand=strand),
                        _SimpleLocation(0, y, strand=strand),
                    )
                )
            else:
                raise err

        sf = _SeqFeature(location, type=type_, qualifiers=qualifiers)

        self.features.append(sf)

        """
        In [11]: a.seq.translate()
        Out[11]: Seq('K', ExtendedIUPACProtein())
        """

    def list_features(self):
        """Print ASCII table with all features.

        Examples
        --------
        >>> from pydna.seq import Seq
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord(Seq("atgtaa"))
        >>> a.add_feature(2,4)
        >>> print(a.list_features())
        +-----+---------------+-----+-----+-----+-----+------+------+
        | Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |
        +-----+---------------+-----+-----+-----+-----+------+------+
        |   0 | L:ft2         | --> | 2   | 4   |   2 | misc |  no  |
        +-----+---------------+-----+-----+-----+-----+------+------+
        """
        x = _PrettyTable(["Ft#", "Label or Note", "Dir", "Sta", "End", "Len", "type", "orf?"])
        x.align["Ft#"] = "r"  # Left align
        x.align["Label or Note"] = "l"  # Left align
        x.align["Len"] = "r"
        x.align["Sta"] = "l"
        x.align["End"] = "l"
        x.align["type"] = "l"
        x.padding_width = 1  # One space between column edges and contents
        for i, sf in enumerate(self.features):
            try:
                lbl = sf.qualifiers["label"]
            except KeyError:
                try:
                    lbl = sf.qualifiers["note"]
                except KeyError:
                    lbl = "nd"
                else:
                    lbl = "N:{}".format(" ".join(lbl).strip())
            else:
                lbl = "L:{}".format(" ".join(lbl).strip())
            x.add_row(
                [
                    i,
                    lbl[:16],
                    {1: "-->", -1: "<--", 0: "---", None: "---"}[sf.location.strand],
                    sf.location.start,
                    sf.location.end,
                    len(sf),
                    sf.type,
                    {True: "yes", False: "no"}[
                        self.extract_feature(i).isorf() or self.extract_feature(i).reverse_complement().isorf()
                    ],
                ]
            )
        return x

    def extract_feature(self, n):
        """Extract feature and return a new SeqRecord object.

        Parameters
        ----------
        n  : int
        Indicates the feature to extract

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> b=a.extract_feature(0)
        >>> b
        SeqRecord(seq=Seq('gt'), id='ft2', name='part_name',
                  description='description', dbxrefs=[])
        """
        return self.features[n].extract(self)

    def sorted_features(self):
        """Return a list of the features sorted by start position.

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.add_feature(3,4)
        >>> a.add_feature(2,4)
        >>> print(a.features)
        [SeqFeature(SimpleLocation(ExactPosition(3), ExactPosition(4),
                                   strand=1),
                    type='misc', qualifiers=...),
         SeqFeature(SimpleLocation(ExactPosition(2), ExactPosition(4),
                                   strand=1),
                    type='misc', qualifiers=...)]
        >>> print(a.sorted_features())
        [SeqFeature(SimpleLocation(ExactPosition(2), ExactPosition(4),
                                   strand=1),
                    type='misc', qualifiers=...),
         SeqFeature(SimpleLocation(ExactPosition(3), ExactPosition(4),
                                   strand=1),
                    type='misc', qualifiers=...)]
        """
        return sorted(self.features, key=lambda x: x.location.start)

    def seguid(self):
        """Return the url safe SEGUID [#]_ for the sequence.

        This checksum is the same as seguid but with base64.urlsafe
        encoding instead of the normal base 64. This means that
        the characters + and / are replaced with - and _ so that
        the checksum can be a part of and URL or a filename.

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("gattaca")
        >>> a.seguid() # original seguid is +bKGnebMkia5kNg/gF7IORXMnIU
        'lsseguid=tp2jzeCM2e3W4yxtrrx09CMKa_8'

        References
        ----------
        .. [#] http://wiki.christophchamp.com/index.php/SEGUID
        """
        return self.seq.seguid()

    def comment(self, newcomment=""):
        """docstring."""
        result = self.annotations.get("comment", "")
        if newcomment:
            self.annotations["comment"] = (result + "\n" + newcomment).strip()
            result = _pretty_str(self.annotations["comment"])
        return result

    def datefunction():
        """docstring."""
        return datetime.datetime.now().replace(microsecond=0).isoformat()

    def stamp(self, now=datefunction, tool="pydna", separator=" ", comment=""):
        """Add seguid checksum to COMMENTS sections

        The checksum is stored in object.annotations["comment"].
        This shows in the COMMENTS section of a formatted genbank file.

        For blunt linear sequences:

        ``SEGUID <seguid>``

        For circular sequences:

        ``cSEGUID <seguid>``

        Fore linear sequences which are not blunt:

        ``lSEGUID <seguid>``


        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a = SeqRecord("aa")
        >>> a.stamp()
        'lsseguid=gBw0Jp907Tg_yX3jNgS4qQWttjU'
        >>> a.annotations["comment"][:41]
        'pydna lsseguid=gBw0Jp907Tg_yX3jNgS4qQWttj'
        """
        chksum = self.seq.seguid()
        oldcomment = self.annotations.get("comment", "")
        oldstamp = _re.findall(r"..seguid=\S{27}", oldcomment)
        if oldstamp and oldstamp[0] == chksum:
            return _pretty_str(oldstamp[0])
        elif oldstamp:
            _warn(
                f"Stamp change.\nNew: {chksum}\nOld: {oldstamp[0]}",
                _PydnaWarning,
            )
        self.annotations["comment"] = (f"{oldcomment}\n" f"{tool} {chksum} {now()} {comment}").strip()
        return _pretty_str(chksum)

    def lcs(self, other, *args, limit=25, **kwargs):
        """Return the longest common substring between the sequence.

        and another sequence (other). The other sequence can be a string,
        Seq, SeqRecord, Dseq or DseqRecord.
        The method returns a SeqFeature with type "read" as this method
        is mostly used to map sequence reads to the sequence. This can be
        changed by passing a type as keyword with some other string value.

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a = SeqRecord("GGATCC")
        >>> a.lcs("GGATCC", limit=6)
        SeqFeature(SimpleLocation(ExactPosition(0),
                                  ExactPosition(6), strand=1),
                                  type='read',
                                  qualifiers=...)
        >>> a.lcs("GATC", limit=4)
        SeqFeature(SimpleLocation(ExactPosition(1),
                                  ExactPosition(5), strand=1),
                                  type='read',
                                  qualifiers=...)
        >>> a = SeqRecord("CCCCC")
        >>> a.lcs("GGATCC", limit=6)
        SeqFeature(None)

        """
        # longest_common_substring
        # https://biopython.org/wiki/ABI_traces
        if hasattr(other, "seq"):
            r = other.seq
            if hasattr(r, "watson"):
                r = str(r.watson).lower()
            else:
                r = str(r).lower()
        else:
            r = str(other.lower())

        olaps = _common_sub_strings(str(self.seq).lower(), r, limit=limit or 25)

        try:
            start_in_self, start_in_other, length = olaps.pop(0)
        except IndexError:
            result = _SeqFeature()
        else:
            label = "sequence" if not hasattr(other, "name") else other.name
            result = _SeqFeature(
                _SimpleLocation(start_in_self, start_in_self + length, strand=1),
                type=kwargs.get("type") or "read",
                qualifiers={
                    "label": [kwargs.get("label") or label],
                    "ApEinfo_fwdcolor": ["#DAFFCF"],
                    "ApEinfo_revcolor": ["#DFFDFF"],
                },
            )
        return result

    def gc(self):
        """Return GC content."""
        return self.seq.gc()

    def cai(self, organism="sce"):
        """docstring."""
        return self.seq.cai(organism=organism)

    def rarecodons(self, organism="sce"):
        """docstring."""
        sfs = []
        for slc in self.seq.rarecodons(organism):
            cdn = self.seq._data[slc].decode("ASCII")
            sfs.append(
                _SeqFeature(
                    _SimpleLocation(slc.start, slc.stop),
                    type=f"rare_codon_{organism}",
                    qualifiers={"label": [cdn]},
                )
            )
        return sfs

    def startcodon(self, organism="sce"):
        """docstring."""
        return self.seq.startcodon()

    def stopcodon(self, organism="sce"):
        """docstring."""
        return self.seq.stopcodon()

    def express(self, organism="sce"):
        """docstring."""
        return self.seq.express()

    def copy(self):
        """docstring."""
        return _copy(self)

    def __lt__(self, other):
        """docstring."""
        try:
            return str(self.seq) < str(other.seq)
        except AttributeError:
            # I don't know how to compare to other
            return NotImplemented

    def __gt__(self, other):
        """docstring."""
        try:
            return str(self.seq) > str(other.seq)
        except AttributeError:
            # I don't know how to compare to other
            return NotImplemented

    def __eq__(self, other):
        """docstring."""
        try:
            if self.seq == other.seq and str(self.__dict__) == str(other.__dict__):
                return True
        except AttributeError:
            pass
        return False

    def __ne__(self, other):
        """docstring."""
        return not self.__eq__(other)

    def __hash__(self):
        """__hash__ must be based on __eq__."""
        return hash((str(self.seq).lower(), str(tuple(sorted(self.__dict__.items())))))

    def __str__(self):
        """docstring."""
        return _pretty_str(super().__str__())

    def __repr__(self):
        """docstring."""
        return _pretty_str(super().__repr__())

    def __format__(self, format):
        """docstring."""

        def removeprefix(text, prefix):
            """Until Python 3.8 is dropped, then use str.removeprefix."""
            if text.startswith(prefix):
                return text[len(prefix) :]
            return text

        if format == "pydnafasta":
            return _pretty_str(
                f">{self.id} {len(self)} bp {dict(((True, 'circular'), (False, 'linear')))[self.seq.circular]}\n{str(self.seq)}\n"
            )
        if format == "primer":
            return _pretty_str(
                f">{self.id} {len(self)}-mer{removeprefix(self.description, self.name).strip()}\n{str(self.seq)}\n"
            )
        return _pretty_str(super().__format__(format))

    def __add__(self, other):
        """docstring."""
        answer = super().__add__(other)
        if answer.name == "<unknown name>":
            answer.name = "name"
        if answer.id == "<unknown id>":
            answer.id = "id"
        if answer.description == "<unknown description>":
            answer.description = "description"
        answer = type(self).from_Bio_SeqRecord(answer)
        return answer

    def __getitem__(self, index):
        """docstring."""
        from pydna.utils import (
            identifier_from_string as _identifier_from_string,
        )  # TODO: clean this up

        answer = super().__getitem__(index)
        if len(answer) < 2:
            return answer
        identifier = "part_{id}".format(id=self.id)
        if answer.features:
            sf = max(answer.features, key=len)  # default
            if "label" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["label"])
            elif "note" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["note"])
        answer.id = _identifier_from_string(identifier)[:16]
        answer.name = _identifier_from_string(f"part_{self.name}")[:16]
        return answer

    def __bool__(self):
        """Boolean value of an instance of this class (True).

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a SeqRecord always evaluated as True.

        Note that in comparison, a Seq object will evaluate to False if it
        has a zero length sequence.

        WARNING: The SeqRecord may in future evaluate to False when its
        sequence is of zero length (in order to better match the Seq
        object behaviour)!
        """
        return bool(self.seq)

    def dump(self, filename, protocol=None):
        """docstring."""
        from pathlib import Path

        pth = Path(filename)
        if not pth.suffix:
            pth = pth.with_suffix(".pickle")
        with open(pth, "wb") as f:
            _pickle.dump(self, f, protocol=protocol)
        return _pretty_str(pth)


class ProteinSeqRecord(SeqRecord):

    def reverse_complement(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    rc = reverse_complement

    def isorf(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def gc(self):
        raise NotImplementedError("Not defined for protein.")

    def cai(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def rarecodons(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def startcodon(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def stopcodon(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def express(self, *args, **kwargs):
        raise NotImplementedError("Not defined for protein.")

    def __format__(self, format):
        """docstring."""
        return _pretty_str(_SeqRecord.__format__(self, format))


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=(doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE))
    _os.environ["pydna_cached_funcs"] = cached
