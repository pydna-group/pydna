#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provides the :class:`Dseqrecord` class, for handling double stranded
DNA sequences. The Dseqrecord holds sequence information in the form of a :class:`pydna.dseq.Dseq`
object. The Dseq and Dseqrecord classes are subclasses of Biopythons
Seq and SeqRecord classes, respectively.

The Dseq and Dseqrecord classes support the notion of circular and linear DNA topology.
"""
from Bio.Restriction import RestrictionBatch as _RestrictionBatch
from Bio.Restriction import CommOnly
from pydna.dseq import Dseq as _Dseq
from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import flatten as _flatten, location_boundaries as _location_boundaries

# from pydna.utils import memorize as _memorize
from pydna.utils import rc as _rc
from pydna.utils import shift_location as _shift_location
from pydna.utils import shift_feature as _shift_feature
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from Bio.SeqFeature import SeqFeature as _SeqFeature
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from pydna.seqrecord import SeqRecord as _SeqRecord
from Bio.Seq import translate as _translate
from pydna.utils import identifier_from_string as _identifier_from_string
import copy as _copy
import operator as _operator
import os as _os
import re as _re
import time as _time
import datetime as _datetime


import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


try:
    from IPython.display import display_html as _display_html
except ImportError:

    def _display_html(item, raw=None):
        return item


class Dseqrecord(_SeqRecord):
    """Dseqrecord is a double stranded version of the Biopython SeqRecord [#]_ class.
    The Dseqrecord object holds a Dseq object describing the sequence.
    Additionally, Dseqrecord hold meta information about the sequence in the
    from of a list of SeqFeatures, in the same way as the SeqRecord does.

    The Dseqrecord can be initialized with a string, Seq, Dseq, SeqRecord
    or another Dseqrecord. The sequence information will be stored in a
    Dseq object in all cases.

    Dseqrecord objects can be read or parsed from sequences in FASTA, EMBL or Genbank formats.
    See the :mod:`pydna.readers` and :mod:`pydna.parsers` modules for further information.

    There is a short representation associated with the Dseqrecord.
    ``Dseqrecord(-3)`` represents a linear sequence of length 2
    while ``Dseqrecord(o7)``
    represents a circular sequence of length 7.

    Dseqrecord and Dseq share the same concept of length. This length can be larger
    than each strand alone if they are staggered as in the example below.

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG




    Parameters
    ----------
    record  : string, Seq, SeqRecord, Dseq or other Dseqrecord object
        This data will be used to form the seq property

    circular : bool, optional
        True or False reflecting the shape of the DNA molecule

    linear : bool, optional
        True or False reflecting the shape of the DNA molecule


    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> a=Dseqrecord("aaa")
    >>> a
    Dseqrecord(-3)
    >>> a.seq
    Dseq(-3)
    aaa
    ttt
    >>> from pydna.seq import Seq
    >>> b=Dseqrecord(Seq("aaa"))
    >>> b
    Dseqrecord(-3)
    >>> b.seq
    Dseq(-3)
    aaa
    ttt
    >>> from Bio.SeqRecord import SeqRecord
    >>> c=Dseqrecord(SeqRecord(Seq("aaa")))
    >>> c
    Dseqrecord(-3)
    >>> c.seq
    Dseq(-3)
    aaa
    ttt

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    """

    def __init__(
        self,
        record,
        *args,
        circular=None,
        n=5e-14,  # mol ( = 0.05 pmol)
        **kwargs,
    ):
        _module_logger.info("### Dseqrecord initialized ###")
        _module_logger.info("argument circular = %s", circular)
        _module_logger.info("circular = %s", circular)

        if isinstance(record, str):
            _module_logger.info("record is a string")
            super().__init__(
                _Dseq.from_string(
                    record,
                    # linear=linear,
                    circular=bool(circular),
                ),
                *args,
                **kwargs,
            )

        # record is a Dseq object ?
        elif hasattr(record, "watson"):
            if circular is False:
                record = record[:]
            elif circular is True:
                record = record.looped()
            _module_logger.info("record is a Dseq object")
            super().__init__(record, *args, **kwargs)

        # record is a Bio.Seq object ?
        elif hasattr(record, "transcribe"):
            _module_logger.info("record is a Seq object")
            super().__init__(
                _Dseq(
                    str(record),
                    # linear=linear,
                    circular=bool(circular),
                ),
                *args,
                **kwargs,
            )

        # record is a Bio.SeqRecord or Dseqrecord object ?
        elif hasattr(record, "features"):
            _module_logger.info("record is a Bio.SeqRecord or Dseqrecord object")
            for key, value in list(record.__dict__.items()):
                setattr(self, key, value)
            self.letter_annotations = {}
            # record.seq is a Dseq object ?
            if hasattr(record.seq, "watson"):
                new_seq = _copy.copy(record.seq)
                if circular is False:
                    new_seq = new_seq[:]
                elif circular is True:
                    new_seq = new_seq.looped()
                self.seq = new_seq
            # record.seq is Bio.SeqRecord object ?
            else:
                self.seq = _Dseq(
                    str(record.seq),
                    # linear=linear,
                    circular=bool(circular),
                )
        else:
            raise ValueError("don't know what to do with {}".format(record))

        self.map_target = None
        self.n = n  # amount, set to 5E-14 which is 5 pmols
        self.annotations.update({"molecule_type": "DNA"})

    @classmethod
    def from_string(
        cls,
        record: str = "",
        *args,
        # linear=True,
        circular=False,
        n=5e-14,
        **kwargs,
    ):
        """docstring."""
        # def from_string(cls, record:str="", *args,
        # linear=True, circular=False, n = 5E-14, **kwargs):
        obj = cls.__new__(cls)  # Does not call __init__
        obj._per_letter_annotations = {}
        obj.seq = _Dseq.quick(
            record,
            _rc(record),
            ovhg=0,
            # linear=linear,
            circular=circular,
        )
        obj.id = _pretty_str("id")
        obj.name = _pretty_str("name")
        obj.description = _pretty_str("description")
        obj.dbxrefs = []
        obj.annotations = {"molecule_type": "DNA"}
        obj.features = []
        obj.map_target = None
        obj.n = n
        obj.__dict__.update(kwargs)
        return obj

    @classmethod
    def from_SeqRecord(
        cls,
        record: _SeqRecord,
        *args,
        circular=None,
        n=5e-14,
        **kwargs,
    ):
        obj = cls.__new__(cls)  # Does not call __init__
        obj._per_letter_annotations = record._per_letter_annotations
        obj.id = record.id
        obj.name = record.name
        obj.description = record.description
        obj.dbxrefs = record.dbxrefs
        obj.annotations = {"molecule_type": "DNA"}
        obj.annotations.update(record.annotations)
        obj.features = record.features
        obj.map_target = None
        obj.n = n
        if circular is None:
            circular = record.annotations.get("topology") == "circular"
        obj.seq = _Dseq.quick(str(record.seq), _rc(str(record.seq)), ovhg=0, circular=circular)
        return obj

    @property
    def circular(self):
        """The circular property can not be set directly.
        Use :meth:`looped`"""
        return self.seq.circular

    def m(self):
        """This method returns the mass of the DNA molecule in grams. This is
        calculated as the product between the molecular weight of the Dseq object
        and the"""
        return self.seq.mw() * self.n  # Da(g/mol) * mol = g

    def extract_feature(self, n):
        """Extracts a feature and creates a new Dseqrecord object.

        Parameters
        ----------
        n  : int
            Indicates the feature to extract

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> b=a.extract_feature(0)
        >>> b
        Dseqrecord(-2)
        >>> b.seq
        Dseq(-2)
        gt
        ca

        """
        return super().extract_feature(n)

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
        [SeqFeature(SimpleLocation(ExactPosition(2), ExactPosition(4), strand=1), type='misc', qualifiers=...)]
        """
        if x and y and self.circular and x > y:
            pass
        else:
            super().add_feature(x, y, seq, type_, strand=strand, *args, **kwargs)
            return

        qualifiers = {}
        qualifiers.update(kwargs)

        location = _CompoundLocation(
            (
                _SimpleLocation(x, self.seq.length, strand=strand),
                _SimpleLocation(0, y, strand=strand),
            )
        )

        sf = _SeqFeature(location, type=type_, qualifiers=qualifiers)

        if "label" not in qualifiers:
            qualifiers["label"] = [f"ft{len(location)}"]

        if sf.extract(self).isorf():
            qualifiers["label"] = [f"orf{len(location)}"]

        self.features.append(sf)

    def seguid(self):
        """Url safe SEGUID for the sequence.

        This checksum is the same as seguid but with base64.urlsafe
        encoding instead of the normal base64. This means that
        the characters + and / are replaced with - and _ so that
        the checksum can be part of a URL.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a = Dseqrecord("aa")
        >>> a.seguid()
        'ldseguid=TEwydy0ugvGXh3VJnVwgtxoyDQA'

        """
        return self.seq.seguid()

    def looped(self):
        """
        Circular version of the Dseqrecord object.

        The underlying linear Dseq object has to have compatible ends.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaa")
        >>> a
        Dseqrecord(-3)
        >>> b=a.looped()
        >>> b
        Dseqrecord(o3)
        >>>

        See Also
        --------
        pydna.dseq.Dseq.looped
        """
        new = _copy.copy(self)
        # for key, value in list(self.__dict__.items()):
        #     setattr(new, key, value)
        new._seq = self.seq.looped()
        five_prime = self.seq.five_prime_end()
        for fn, fo in zip(new.features, self.features):
            if five_prime[0] == "5'":
                pass
                # fn.location = fn.location + self.seq.ovhg
            elif five_prime[0] == "3'":
                fn.location = fn.location + (-self.seq.ovhg)
            if fn.location.start < 0:
                loc1 = _SimpleLocation(len(new) + fn.location.start, len(new), strand=fn.location.strand)
                loc2 = _SimpleLocation(0, fn.location.end, strand=fn.location.strand)
                fn.location = _CompoundLocation([loc1, loc2])

            if fn.location.end > len(new):
                loc1 = _SimpleLocation(fn.location.start, len(new), strand=fn.location.strand)
                loc2 = _SimpleLocation(0, fn.location.end - len(new), strand=fn.location.strand)
                fn.location = _CompoundLocation([loc1, loc2])

            fn.qualifiers = fo.qualifiers

        return new

    def tolinear(self):  # pragma: no cover
        """
        Returns a linear, blunt copy of a circular Dseqrecord object. The
        underlying Dseq object has to be circular.

        This method is deprecated, use slicing instead. See example below.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaa", circular = True)
        >>> a
        Dseqrecord(o3)
        >>> b=a[:]
        >>> b
        Dseqrecord(-3)
        >>>

        """
        import warnings as _warnings
        from pydna import _PydnaDeprecationWarning

        _warnings.warn(
            "tolinear method is obsolete; " "please use obj[:] " "instead of obj.tolinear().",
            _PydnaDeprecationWarning,
        )
        new = _copy.copy(self)
        for key, value in list(self.__dict__.items()):
            setattr(new, key, value)
        # new._seq = self.seq.tolinear()
        for fn, fo in zip(new.features, self.features):
            fn.qualifiers = fo.qualifiers

        return new

    def terminal_transferase(self, nucleotides="a"):
        """docstring."""
        newseq = _copy.deepcopy(self)
        newseq.seq = self.seq.terminal_transferase(nucleotides)
        for feature in newseq.features:
            feature.location += len(nucleotides)
        return newseq

    def format(self, f="gb"):
        """Returns the sequence as a string using a format supported by Biopython
        SeqIO [#]_. Default is "gb" which is short for Genbank.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> x=Dseqrecord("aaa")
        >>> x.annotations['date'] = '02-FEB-2013'
        >>> x
        Dseqrecord(-3)
        >>> print(x.format("gb"))
        LOCUS       name                       3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  description.
        ACCESSION   id
        VERSION     id
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //


        References
        ----------

        .. [#] http://biopython.org/wiki/SeqIO


        """

        record = _copy.deepcopy(self)
        if f in ("genbank", "gb") and self.circular:
            record.annotations["topology"] = "circular"
        else:
            record.annotations["topology"] = "linear"
        return _SeqRecord.format(record, f).strip()

    def write(self, filename=None, f="gb"):
        """Writes the Dseqrecord to a file using the format f, which must
        be a format supported by Biopython SeqIO for writing [#]_. Default
        is "gb" which is short for Genbank. Note that Biopython SeqIO reads
        more formats than it writes.

        Filename is the path to the file where the sequece is to be
        written. The filename is optional, if it is not given, the
        description property (string) is used together with the format.

        If obj is the Dseqrecord object, the default file name will be:

        ``<obj.locus>.<f>``

        Where <f> is "gb" by default. If the filename already exists and
        AND the sequence it contains is different, a new file name will be
        used so that the old file is not lost:

        ``<obj.locus>_NEW.<f>``

        References
        ----------
        .. [#] http://biopython.org/wiki/SeqIO

        """
        msg = ""
        if not filename:
            filename = "{name}.{type}".format(name=self.locus, type=f)
            # generate a name if no name was given
        # if not isinstance(filename, str):  # is filename a string???
        #     raise ValueError("filename has to be a string, got", type(filename))
        name, ext = _os.path.splitext(filename)
        msg = f"<font face=monospace><a href='{filename}' target='_blank'>{filename}</a></font><br>"
        if not _os.path.isfile(filename):
            with open(filename, "w", encoding="utf8") as fp:
                fp.write(self.format(f))
        else:
            from pydna.readers import read

            old_file = read(filename)

            if self.seq != old_file.seq:
                # If new sequence is different, the old file is
                # renamed with "_OLD_" suffix:
                oldmtime = _datetime.datetime.fromtimestamp(_os.path.getmtime(filename)).isoformat()
                tstmp = int(_time.time() * 1_000_000)
                old_filename = f"{name}_OLD_{tstmp}{ext}"
                _os.rename(filename, old_filename)
                with open(filename, "w", encoding="utf8") as fp:
                    fp.write(self.format(f))
                newmtime = _datetime.datetime.fromtimestamp(_os.path.getmtime(filename)).isoformat()
                msg = f"""
                <table style="padding:10px 10px;
                word-break:normal;
                border-color:#fe0000;
                border-collapse:collapse;
                border-spacing:1;
                font-family:monospace;
                font-size:large;
                font-weight:bold;
                text-align:left;
                border: 5px solid red;">
                <thead>
                  <tr style="color:#0000FF;border: 1px solid;text-align:left;">
                    <th style="color:#fe0000;border: 1px solid;text-align:center;font-size:xxx-large;text-align:left;">&#9888</th>
                    <th style="color:#f56b00;border: 1px solid;text-align:left;" colspan="2">Sequence change</th>
                  </tr>
                </thead>
                <tbody>
                  <tr style="color:#0000FF;border: 1px solid;text-align:left;">
                    <td>Filename</td>
                    <td style="color:#fe0000;border: 1px solid;text-align:left;"><a href='{filename}' target='_blank'>{filename}</a></td>
                    <td style="color:#32cb00;border: 1px solid;text-align:left;"><a href='{old_filename}' target='_blank'>{old_filename}</a></td>
                  </tr>
                  <tr style="color:#0000FF;border: 1px solid;text-align:left;">
                    <td >Saved</td>
                    <td style="color:#fe0000;border: 1px solid;text-align:left;">{newmtime}</td>
                    <td style="color:#32cb00;border: 1px solid;text-align:left;">{oldmtime}</td>
                  </tr>
                  <tr style="color:#0000FF;border: 1px solid;text-align:left;">
                    <td>Length</td>
                    <td style="color:#fe0000;border: 1px solid;text-align:left;">{len(self)}</td>
                    <td style="color:#32cb00;border: 1px solid;text-align:left;">{len(old_file)}</td>
                  </tr>
                  <tr style="color:#0000FF;border: 1px solid;text-align:left;">
                    <td>SEGUID</td>
                    <td style="color:#fe0000;border: 1px solid;text-align:left;">{self.seguid()}</td>
                    <td style="color:#32cb00;border: 1px solid;text-align:left;">{old_file.seguid()}</td>
                  </tr>
                </tbody>
                </table>
                """
            elif "seguid" in old_file.annotations.get("comment", ""):
                pattern = r"(ldseguid|cdseguid)-(\S{27})(_[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{6}){0,1}"
                # seguid=NNNNNNNNNNNNNNNNNNNNNNNNNNN_2020-10-10T11:11:11.111111
                oldstamp = _re.search(pattern, old_file.description)
                newstamp = _re.search(pattern, self.description)
                newdescription = self.description
                if oldstamp and newstamp:
                    if oldstamp.group(0)[:35] == newstamp.group(0)[:35]:
                        newdescription = newdescription.replace(newstamp.group(0), oldstamp.group(0))
                elif oldstamp:
                    newdescription += " " + oldstamp.group(0)
                newobj = _copy.copy(self)
                newobj.description = newdescription

                with open(filename, "w", encoding="utf8") as fp:
                    fp.write(newobj.format(f))
            else:
                with open(filename, "w", encoding="utf8") as fp:
                    fp.write(self.format(f))
        return _display_html(msg, raw=True)

    def find(self, other):
        # TODO allow strings, seqs, seqrecords or Dseqrecords
        # TODO check for linearity of other, raise exception if not
        # TODO add tests and docstring for this method
        o = str(other.seq).upper()

        if not self.circular:
            s = str(self.seq).upper()
        else:
            # allow wrapping around origin
            s = str(self.seq).upper() + str(self.seq).upper()[: len(other) - 1]
        return s.find(o)

    def __str__(self):
        return ("Dseqrecord\n" "circular: {}\n" "size: {}\n").format(self.circular, len(self)) + _SeqRecord.__str__(
            self
        )

    def __contains__(self, other):
        if other.lower() in str(self.seq).lower():
            return True
        else:
            s = self.seq.watson.replace(" ", "")
            ln = len(s)
            spc = 3 - ln % 3 if ln % 3 else 0
            s = "n" * spc + s + "nnn"
            for frame in range(3):
                if other.lower() in _translate(s[frame : frame + spc + ln]).lower():
                    return True
        return False

    def find_aminoacids(self, other):
        """
        >>> from pydna.dseqrecord import Dseqrecord
        >>> s=Dseqrecord("atgtacgatcgtatgctggttatattttag")
        >>> s.seq.translate()
        Seq('MYDRMLVIF*')
        >>> "RML" in s
        True
        >>> "MMM" in s
        False
        >>> s.seq.rc().translate()
        Seq('LKYNQHTIVH')
        >>> "QHT" in s.rc()
        True
        >>> "QHT" in s
        False
        >>> slc = s.find_aa("RML")
        >>> slc
        slice(9, 18, None)
        >>> s[slc]
        Dseqrecord(-9)
        >>> code = s[slc].seq
        >>> code
        Dseq(-9)
        cgtatgctg
        gcatacgac
        >>> code.translate()
        Seq('RML')
        """
        other = str(other).lower()
        assert self.seq.watson == "".join(self.seq.watson.split())
        s = self.seq.watson
        ln = len(s)
        spc = 3 - ln % 3 if ln % 3 else 0
        s = s + "n" * spc + "nnn"
        start = None
        for frame in range(3):
            try:
                start = _translate(s[frame : frame + ln + spc]).lower().index(other)
                break
            except ValueError:
                pass
        oh = self.seq.ovhg if self.seq.ovhg > 0 else 0
        if start is None:
            return None  # TODO return an emoty slice or False...?
        else:
            return slice(frame + start * 3 + oh, frame + (start + len(other)) * 3 + oh)

    find_aa = find_aminoacids

    def map_trace_files(self, pth, limit=25):  # TODO allow path-like objects
        import glob

        traces = []
        for name in glob.glob(pth):
            trace = SeqIO.read(name, "abi").lower()
            trace.annotations["filename"] = trace.fname = name
            traces.append(trace)
        if not traces:
            raise ValueError("No trace files found in {}".format(pth))
        if hasattr(self.map_target, "step"):
            area = self.map_target
        elif hasattr(self.map_target, "extract"):
            area = slice(self.map_target.location.start, self.map_target.location.end)
        else:
            area = None  # TODO allow other objects as well and do some checks on map target

        if area:
            self.matching_reads = []
            self.not_matching_reads = []
            target = str(self[area].seq).lower()
            target_rc = str(self[area].seq.rc()).lower()
            for trace in traces:
                if target in str(trace.seq) or target_rc in str(trace.seq):
                    self.matching_reads.append(trace)
                else:
                    self.not_matching_reads.append(trace)
            reads = self.matching_reads
        else:
            self.matching_reads = None
            self.not_matching_reads = None
            reads = traces

        matching_reads = []

        for read_ in reads:
            matches = _common_sub_strings(str(self.seq).lower(), str(read_.seq), limit)

            if not matches:
                continue

            if len(matches) > 1:
                newmatches = [
                    matches[0],
                ]
                for i, x in enumerate(matches[1:]):
                    g, f, h = matches[i]
                    if g + h < x[0] and f + h < x[1]:
                        newmatches.append(x)
            else:  # len(matches)==1
                newmatches = matches

            matching_reads.append(read_)

            if len(newmatches) > 1:
                ms = []
                for m in newmatches:
                    ms.append(_SimpleLocation(m[0], m[0] + m[2]))
                loc = _CompoundLocation(ms)
            else:
                a, b, c = newmatches[0]
                loc = _SimpleLocation(a, a + c)

            self.features.append(
                _SeqFeature(
                    loc,
                    qualifiers={"label": [read_.annotations["filename"]]},
                    type="trace",
                )
            )

        return [x.annotations["filename"] for x in matching_reads]

    def __repr__(self):
        return "Dseqrecord({}{})".format({True: "-", False: "o"}[not self.circular], len(self))

    def _repr_pretty_(self, p, cycle):
        p.text("Dseqrecord({}{})".format({True: "-", False: "o"}[not self.circular], len(self)))

    def __add__(self, other):
        if hasattr(other, "seq") and hasattr(other.seq, "watson"):
            other = _copy.deepcopy(other)
            other_five_prime = other.seq.five_prime_end()
            if other_five_prime[0] == "5'":
                # add other.seq.ovhg
                for f in other.features:
                    f.location = f.location + other.seq.ovhg
            elif other_five_prime[0] == "3'":
                # subtract other.seq.ovhg (sign change)
                for f in other.features:
                    f.location = f.location + (-other.seq.ovhg)

            answer = Dseqrecord(_SeqRecord.__add__(self, other))
            answer.n = min(self.n, other.n)
        else:
            answer = Dseqrecord(_SeqRecord.__add__(self, Dseqrecord(other)))
            answer.n = self.n
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseqrecord by non-int of type {}".format(type(number)))
        if self.circular:
            raise TypeError("TypeError: can't multiply circular Dseqrecord.")
        if number > 0:
            new = _copy.deepcopy(self)
            for i in range(1, number):
                new += self
            new._per_letter_annotations = self._per_letter_annotations
            return new
        else:
            return self.__class__("")

    def __getitem__(self, sl):
        """docstring."""
        answer = Dseqrecord(_copy.copy(self))
        answer.seq = self.seq.__getitem__(sl)
        # answer.seq.alphabet = self.seq.alphabet
        # breakpoint()
        sl_start = sl.start or 0  # 6
        sl_stop = sl.stop or len(self.seq)  # 1

        if not self.circular or sl_start < sl_stop:
            # TODO: special case for sl_end == 0 in circular sequences
            # related to https://github.com/pydna-group/pydna/issues/161
            if self.circular and sl.stop == 0:
                sl = slice(sl.start, len(self.seq), sl.step)
            answer.features = super().__getitem__(sl).features
        elif self.circular and sl_start > sl_stop:
            answer.features = self.shifted(sl_start).features
            # origin-spanning features should only be included after shifting
            # in cases where the slice comprises the entire sequence, but then
            # sl_start == sl_stop and the second condition is not met
            answer.features = [
                f
                for f in answer.features
                if (
                    _location_boundaries(f.location)[1] <= answer.seq.length
                    and _location_boundaries(f.location)[0] < _location_boundaries(f.location)[1]
                )
            ]

        elif self.circular and sl_start == sl_stop:
            cut = ((sl_start, 0), None)
            return self.apply_cut(cut, cut)
        else:
            answer = Dseqrecord("")
        identifier = "part_{id}".format(id=self.id)
        if answer.features:
            sf = max(answer.features, key=len)  # default
            if "label" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["label"])
            elif "note" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["note"])
        answer.id = _identifier_from_string(identifier)[:16]
        answer.name = _identifier_from_string("part_{name}".format(name=self.name))[:16]
        return answer

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

    def linearize(self, *enzymes):
        """Similar to `:func:cut`.

        Throws an exception if there is not excactly one cut
        i.e. none or more than one digestion products.
        """
        if not self.seq.circular:
            raise TypeError("Can only linearize circular molecules!")
        fragments = self.cut(*enzymes)
        if len(fragments) > 1:
            raise TypeError("More than one fragment is formed!")
        elif not fragments:
            raise TypeError("The enzyme(s) do not cut!")
        answer = fragments[0]
        answer.id = "{name}_lin".format(name=self.name)
        answer.name = answer.id[:16]
        return fragments[0]

    def no_cutters(self, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.no_cutters(batch=batch or CommOnly)

    def unique_cutters(self, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.unique_cutters(batch=batch or CommOnly)

    def once_cutters(self, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.once_cutters(batch=batch or CommOnly)

    def twice_cutters(self, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.twice_cutters(batch=batch or CommOnly)

    def n_cutters(self, n=3, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.n_cutters(n=n, batch=batch or CommOnly)

    def cutters(self, batch: _RestrictionBatch = None):
        """docstring."""
        return self.seq.cutters(batch=batch or CommOnly)

    def number_of_cuts(self, *enzymes):
        """The number of cuts by digestion with the Restriction enzymes
        contained in the iterable."""
        return sum([len(enzyme.search(self.seq)) for enzyme in _flatten(enzymes)])

    def cas9(self, RNA: str):
        """docstring."""
        fragments = []
        result = []
        for target in (self.seq, self.seq.rc()):
            fragments = [self[sl.start : sl.stop] for sl in target.cas9(RNA)]
            result.append(fragments)
        return result

    def reverse_complement(self):
        """Reverse complement.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("ggaatt")
        >>> a
        Dseqrecord(-6)
        >>> a.seq
        Dseq(-6)
        ggaatt
        ccttaa
        >>> a.reverse_complement().seq
        Dseq(-6)
        aattcc
        ttaagg
        >>>

        See Also
        --------
        pydna.dseq.Dseq.reverse_complement

        """
        answer = type(self)(super().reverse_complement())
        answer.name = "{}_rc".format(self.name[:13])
        answer.description = self.description + "_rc"
        answer.id = self.id + "_rc"
        answer.seq.circular = self.seq.circular
        # answer.seq._linear = self.seq.linear
        return answer

    rc = reverse_complement

    # @_memorize("pydna.dseqrecord.Dseqrecord.synced")
    def synced(self, ref, limit=25):
        """This method returns a new circular sequence (Dseqrecord object), which has been rotated
        in such a way that there is maximum overlap between the sequence and
        ref, which may be a string, Biopython Seq, SeqRecord object or
        another Dseqrecord object.

        The reason for using this could be to rotate a new recombinant plasmid so
        that it starts at the same position after cloning. See the example below:


        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("gaat", circular=True)
        >>> a.seq
        Dseq(o4)
        gaat
        ctta
        >>> d = a[2:] + a[:2]
        >>> d.seq
        Dseq(-4)
        atga
        tact
        >>> insert=Dseqrecord("CCC")
        >>> recombinant = (d+insert).looped()
        >>> recombinant.seq
        Dseq(o7)
        atgaCCC
        tactGGG
        >>> recombinant.synced(a).seq
        Dseq(o7)
        gaCCCat
        ctGGGta

        """

        if not self.circular:
            raise TypeError("Only circular DNA can be synced!")

        newseq = _copy.copy(self)

        s = str(self.seq.watson).lower()
        s_rc = str(self.seq.crick).lower()

        if hasattr(ref, "seq"):
            r = ref.seq
            if hasattr(r, "watson"):
                r = str(r.watson).lower()
            else:
                r = str(r).lower()
        else:
            r = str(ref.lower())

        lim = min(limit, limit * (len(s) // limit) + 1)

        c = _common_sub_strings(s + s, r, limit=lim)
        d = _common_sub_strings(s_rc + s_rc, r, limit=lim)

        c = [(x[0], x[2]) for x in c if x[1] == 0]
        d = [(x[0], x[2]) for x in d if x[1] == 0]

        if not c and not d:
            raise TypeError("There is no overlap between sequences!")

        if c:
            start, length = c.pop(0)
        else:
            start, length = 0, 0

        if d:
            start_rc, length_rc = d.pop(0)
        else:
            start_rc, length_rc = 0, 0

        if length_rc > length:
            start = start_rc
            newseq = newseq.rc()

        if start == 0:
            result = newseq
        else:
            result = newseq.shifted(start)
        _module_logger.info("synced")
        return result

    def upper(self):
        """Returns an uppercase copy.
        >>> from pydna.dseqrecord import Dseqrecord
        >>> my_seq = Dseqrecord("aAa")
        >>> my_seq.seq
        Dseq(-3)
        aAa
        tTt
        >>> upper = my_seq.upper()
        >>> upper.seq
        Dseq(-3)
        AAA
        TTT
        >>>


        Returns
        -------
        Dseqrecord
            Dseqrecord object in uppercase


        See also
        --------
        pydna.dseqrecord.Dseqrecord.lower"""

        upper = _copy.deepcopy(self)
        # This is because the @seq.setter methods otherwise sets the _per_letter_annotations to an empty dict
        prev_per_letter_annotation = upper._per_letter_annotations
        upper.seq = upper.seq.upper()
        upper._per_letter_annotations = prev_per_letter_annotation
        return upper

    def lower(self):
        """>>> from pydna.dseqrecord import Dseqrecord
        >>> my_seq = Dseqrecord("aAa")
        >>> my_seq.seq
        Dseq(-3)
        aAa
        tTt
        >>> upper = my_seq.upper()
        >>> upper.seq
        Dseq(-3)
        AAA
        TTT
        >>> lower = my_seq.lower()
        >>> lower
        Dseqrecord(-3)
        >>>

        Returns
        -------
        Dseqrecord
            Dseqrecord object in lowercase

        See also
        --------
        pydna.dseqrecord.Dseqrecord.upper

        """
        lower = _copy.deepcopy(self)
        prev_per_letter_annotation = lower._per_letter_annotations
        lower.seq = lower.seq.lower()
        lower._per_letter_annotations = prev_per_letter_annotation
        return lower

    def orfs(self, minsize=300):
        """docstring."""
        return tuple(Dseqrecord(self[x:y]) for x, y in self.seq.orfs(minsize=minsize))

    def orfs_to_features(self, minsize=300):
        """docstring."""
        features = []
        for strand, s in ((1, self.seq), (-1, self.seq.rc())):
            for x, y in s.orfs(minsize=minsize):
                orf = self[x:y]
                prt = orf.translate()
                features.append(
                    _SeqFeature(
                        _SimpleLocation(x, y, strand=strand),
                        type="CDS",
                        qualifiers={
                            "note": f"{y - x}bp {(y - x) // 3}aa",
                            "checksum": [orf.seguid() + " (DNA)", prt.seguid() + " (protein)"],
                            "codon_start": 1,
                            "transl_table": 11,
                            "translation": str(prt.seq),
                        },
                    )
                )
        return features

    def _copy_to_clipboard(self, sequence_format):
        """docstring."""
        from pyperclip import copy

        copy(self.format(sequence_format))
        return None

    def copy_gb_to_clipboard(self):
        """docstring."""
        self._copy_to_clipboard("gb")
        return None

    def copy_fasta_to_clipboard(self):
        """docstring."""
        self._copy_to_clipboard("fasta")
        return None

    def figure(self, feature=0, highlight="\x1b[48;5;11m", plain="\x1b[0m"):
        """docstring."""
        if self.features:
            f = self.features[feature]
            locations = sorted(self.features[feature].location.parts, key=_SimpleLocation.start.fget)
            strand = f.location.strand
        else:
            locations = [_SimpleLocation(0, 0, 1)]
            strand = 1

        ovhg = self.seq.ovhg + len(self.seq.watson) - len(self.seq.crick)

        w = f"{self.seq.ovhg * chr(32)}{self.seq.watson}{-ovhg * chr(32)}"
        c = f"{-self.seq.ovhg * chr(32)}{self.seq.crick[::-1]}{ovhg * chr(32)}"

        if strand == 1:
            s1, s2 = w, c
        else:
            s1, s2 = c, w

        wfe = [f"{highlight}{s1[part.start:part.end]}{plain}" for part in locations]

        wfe.append("")

        wof = [s1[0 : locations[0].start]]
        for f, s in zip(locations, locations[1:]):
            wof.append(s1[f.end : s.start])
        wof.append(s1[locations[-1].end : len(self)])

        topology = {True: "-", False: "o"}[not self.circular]
        result = f"{self.__class__.__name__}({topology}{len(self)})\n"

        s1 = "".join(f + s for f, s in zip(wof, wfe))

        if strand == 1:
            result += f"{s1}\n{s2}"
        else:
            result += f"{s2}\n{s1}"
        return _pretty_str(result)

    def shifted(self, shift):
        """Circular Dseqrecord with a new origin <shift>.

        This only works on circular Dseqrecords. If we consider the following
        circular sequence:


        | ``GAAAT   <-- watson strand``
        | ``CTTTA   <-- crick strand``

        The T and the G on the watson strand are linked together as well
        as the A and the C of the of the crick strand.

        if ``shift`` is 1, this indicates a new origin at position 1:

        |    new origin at the | symbol:
        |
        | ``G|AAAT``
        | ``C|TTTA``

        new sequence:

        | ``AAATG``
        | ``TTTAC``

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaat",circular=True)
        >>> a
        Dseqrecord(o4)
        >>> a.seq
        Dseq(o4)
        aaat
        ttta
        >>> b=a.shifted(1)
        >>> b
        Dseqrecord(o4)
        >>> b.seq
        Dseq(o4)
        aata
        ttat

        """
        if not self.circular:
            raise TypeError("Sequence is linear, origin can only be " "shifted for circular sequences.\n")
        ln = len(self)
        if not shift % ln:
            return _copy.deepcopy(self)  # shift is a multiple of ln or 0
        else:
            shift %= ln  # 0<=shift<=ln
        newseq = (self.seq[shift:] + self.seq[:shift]).looped()
        newfeatures = _copy.deepcopy(self.features)
        for feature in newfeatures:
            feature.location = _shift_location(feature.location, -shift, ln)
        newfeatures.sort(key=_operator.attrgetter("location.start"))
        answer = _copy.deepcopy(self)
        answer.features = newfeatures
        answer.seq = newseq
        return answer

    def cut(self, *enzymes):
        """Digest a Dseqrecord object with one or more restriction enzymes.

        returns a list of linear Dseqrecords. If there are no cuts, an empty
        list is returned.

        See also :func:`Dseq.cut`
        Parameters
        ----------

        enzymes : enzyme object or iterable of such objects
            A Bio.Restriction.XXX restriction object or iterable of such.

        Returns
        -------
        Dseqrecord_frags : list
            list of Dseqrecord objects formed by the digestion

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("ggatcc")
        >>> from Bio.Restriction import BamHI
        >>> a.cut(BamHI)
        (Dseqrecord(-5), Dseqrecord(-5))
        >>> frag1, frag2 = a.cut(BamHI)
        >>> frag1.seq
        Dseq(-5)
        g
        cctag
        >>> frag2.seq
        Dseq(-5)
        gatcc
            g


        """

        cutsites = self.seq.get_cutsites(*enzymes)
        cutsite_pairs = self.seq.get_cutsite_pairs(cutsites)
        return tuple(self.apply_cut(*cs) for cs in cutsite_pairs)

    def apply_cut(self, left_cut, right_cut):
        dseq = self.seq.apply_cut(left_cut, right_cut)
        # TODO: maybe remove depending on https://github.com/pydna-group/pydna/issues/161

        if left_cut == right_cut:
            # Not really a cut, but to handle the general case
            if left_cut is None:
                features = _copy.deepcopy(self.features)
            else:
                # The features that span the origin if shifting with left_cut, but that do not cross
                # the cut site should be included, and if there is a feature within the cut site, it should
                # be duplicated. See https://github.com/pydna-group/pydna/issues/180 for a practical example.
                #
                # Let's say we are going to open a circular plasmid like below (| inidicate cuts, numbers indicate
                # features)
                #
                #    3333|3
                #    1111
                #     000
                # XXXXatg|YYY
                # XXX|tacYYYY
                #     000
                #     2222
                #
                left_watson, left_crick, left_ovhg = self.seq.get_cut_parameters(left_cut, True)
                initial_shift = left_watson if left_ovhg < 0 else left_crick
                features = self.shifted(initial_shift).features
                # for f in features:
                #     print(f.id, f.location, _location_boundaries(f.location))
                # Here, we have done what's shown below (* indicates the origin).
                # The features 0 and 2 have the right location for the final product:
                #
                #    3*3333
                #    1*111
                # XXXX*atgYYY
                # XXXX*tacYYY
                #      000
                #      2222

                features_need_transfer = [
                    f for f in features if (_location_boundaries(f.location)[1] <= abs(left_ovhg))
                ]
                features_need_transfer = [
                    _shift_feature(f, -abs(left_ovhg), len(self)) for f in features_need_transfer
                ]

                #                                           ^                ^^^^^^^^^
                # Now we have shifted the features that end before the cut (0 and 1, but not 3), as if
                # they referred to the below sequence (* indicates the origin):
                #
                #    1111
                #     000
                # XXXXatg*YYY
                # XXXXtac*YYY
                #
                # The features 0 and 1 would have the right location if the final sequence had the same length
                # as the original one. However, the final product is longer because of the overhang.

                features += [_shift_feature(f, abs(left_ovhg), len(dseq)) for f in features_need_transfer]
                #                             ^                ^^^^^^^^^
                # So we shift back by the same amount in the opposite direction, but this time we pass the
                # length of the final product.
                # print(*features, sep='\n')
                # Features like 3 are removed here
                features = [
                    f
                    for f in features
                    if (
                        _location_boundaries(f.location)[1] <= len(dseq)
                        and _location_boundaries(f.location)[0] <= _location_boundaries(f.location)[1]
                    )
                ]
        else:
            left_watson, left_crick, left_ovhg = self.seq.get_cut_parameters(left_cut, True)
            right_watson, right_crick, right_ovhg = self.seq.get_cut_parameters(right_cut, False)

            left_edge = left_crick if left_ovhg > 0 else left_watson
            right_edge = right_watson if right_ovhg > 0 else right_crick
            features = self[left_edge:right_edge].features

        return Dseqrecord(dseq, features=features)


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    # _os.environ["pydna_cache"] = cache
