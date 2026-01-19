#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the pydna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
# doctest: +IGNORE_EXCEPTION_DETAIL
"""This module provides the :func:`gbtext_clean` function which can clean up broken Genbank files enough to
pass the BioPython Genbank parser

Almost all of this code was lifted from BioJSON (https://github.com/levskaya/BioJSON) by Anselm Levskaya.
The original code was not accompanied by any software licence. This parser is based on pyparsing.

There are some modifications to deal with fringe cases.

The parser first produces JSON as an intermediate format which is then formatted back into a
string in Genbank format.

The parser is not complete, so some fields do not survive the roundtrip (see below).
This should not be a difficult fix. The returned result has two properties,
.jseq which is the intermediate JSON produced by the parser and .gbtext
which is the formatted genbank string."""


import re
import pyparsing as pp

GoodLocus = (
    pp.Literal("LOCUS")
    + pp.Word(pp.alphas + pp.nums + "-_()." + "\\").setResultsName("name")
    + pp.Word(pp.nums).setResultsName("size")
    + pp.Suppress(pp.CaselessLiteral("bp"))
    + pp.Word(pp.alphas + "-").setResultsName("seqtype")
    + (pp.CaselessLiteral("linear") | pp.CaselessLiteral("circular")).setResultsName(
        "topology"
    )
    + pp.Optional(pp.Word(pp.alphas), default="   ").setResultsName("divcode")
    + pp.Regex(r"(\d{2})-(\S{3})-(\d{4})").setResultsName("date")
)

# Older versions of ApE don't include a LOCUS name! Need separate def for this case:
BrokenLocus1 = (
    pp.Literal("LOCUS").setResultsName("name")
    + pp.Word(pp.nums).setResultsName("size")
    + pp.Suppress(pp.CaselessLiteral("bp"))
    + pp.Word(pp.alphas + "-").setResultsName("seqtype")
    + (pp.CaselessLiteral("linear") | pp.CaselessLiteral("circular")).setResultsName(
        "topology"
    )
    + pp.Optional(pp.Word(pp.alphas), default="   ").setResultsName("divcode")
    + pp.Regex(r"(\d{2})-(\S{3})-(\d{4})").setResultsName("date")
)

# LOCUS       YEplac181	5741 bp 	DNA	SYN
BrokenLocus2 = (
    pp.Literal("LOCUS")
    + pp.Word(pp.alphas + pp.nums + "-_()." + "\\").setResultsName("name")
    + pp.Word(pp.nums).setResultsName("size")
    + pp.Suppress(pp.CaselessLiteral("bp"))
    + pp.Word(pp.alphas + "-").setResultsName("seqtype")
    + pp.Optional(
        pp.CaselessLiteral("linear") | pp.CaselessLiteral("circular"),
        default="linear",
    ).setResultsName("topology")
    + pp.Optional(pp.Word(pp.alphas), default="   ").setResultsName("divcode")
    + pp.Regex(r"(\d{2})-(\S{3})-(\d{4})").setResultsName("date")
)

BrokenLocus3 = (
    pp.Literal("LOCUS")
    + pp.Word(pp.alphas + pp.nums + "-_()." + "\\").setResultsName("name")
    + pp.Word(pp.nums).setResultsName("size")
    + pp.Suppress(pp.CaselessLiteral("bp"))
    + pp.Word(pp.alphas + "-").setResultsName("seqtype")
    + pp.Optional(
        pp.CaselessLiteral("linear") | pp.CaselessLiteral("circular"),
        default="linear",
    ).setResultsName("topology")
    + pp.Word(pp.alphas).setResultsName("divcode")
    + pp.Optional(
        pp.Regex(r"(\d{2})-(\S{3})-(\d{4})").setResultsName("date"),
        default="19-MAR-1970",
    ).setResultsName("date")
)

LocusEntry = GoodLocus | BrokenLocus1 | BrokenLocus2 | BrokenLocus3

# ===============================================================================
# Generic Entry

# this catches everything but the FEATURES and SEQUENCE entries, really should add parsing code for
# ACCESSION, COMMENTS, REFERENCE, ORGANISM, etc.
# (Though these entries are generally useless when it comes to hacking on DNA)

# All entries in a genbank file headed by an all-caps title with no space between start-of-line and title
CapWord = pp.Word("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
# after titled line, all subsequent lines have to have at least one space in front of them
# this is how we split up the genbank record
SpacedLine = pp.White(min=1) + pp.CharsNotIn("\n") + pp.LineEnd()
# HeaderLine = CapWord + CharsNotIn("\n") + LineEnd()
GenericEntry = pp.Group(
    CapWord + pp.Combine(pp.CharsNotIn("\n") + pp.LineEnd() + pp.ZeroOrMore(SpacedLine))
).setResultsName("generics", listAllMatches=True)


# ===============================================================================
# Definition Entry
# SuppressedSpacedLine =  Suppress(White(min=1)) + CharsNotIn("\n") + LineEnd()
# DefinitionEntry =  Suppress(Literal("DEFINITION")) + Combine(CharsNotIn("\n") + LineEnd() + ZeroOrMore( SuppressedSpacedLine ))

# ===============================================================================
# GenBank Feature Table Parser

# ==== Genbank Location String Parser
#
# a string of slices w. functional modifiers that go at most two levels deep
# single position is just a number i.e. 23423
# slice is N1..N2  w. N1<N2
# i.e.
# 23..88  --> seq[23:89] in python syntax (genbank uses inclusive slicing)
# 234..555
# complement(234..434) --> rc(seq[234:435])
# join(23..343,454..666,777..999) --> seq[23:344]+seq[454:667]+seq[777:1000]
# complement(join(23..343,454..666,777..999))
# join(complement(34..123),complement(333..565))
#
# additionally the slices can have ambiguous locs like <454..999 or 232..>331
# also note the dumb 34.38 fuzzy slice notation
# i.e. <45..900  says the real feature starts "somewhere" before pos 45
#       45.48 says feature somewhere between bases 45->48
# lot of weird annotations best avoided because they connote ~useless knowledge for synthetic design
#
# if you don't know where something is, don't use it or guess and move on

LPAREN = pp.Suppress("(")
RPAREN = pp.Suppress(")")
SEP = pp.Suppress(pp.Literal(".."))

# recognize numbers w. < & > uncertainty specs, then strip the <> chars to make it fixed
gbIndex = pp.Word(pp.nums + "<>").setParseAction(
    lambda s, l_, t: int(t[0].replace("<", "").replace(">", ""))
)
SimpleSlice = pp.Group(gbIndex + SEP + gbIndex) | pp.Group(gbIndex).setParseAction(
    lambda s, l_, t: [[t[0][0], t[0][0]]]
)

# recursive def for nested function syntax:  f( g(), g() )
complexSlice = pp.Forward()
(
    complexSlice
    << (pp.Literal("complement") | pp.Literal("join"))
    + LPAREN
    + (pp.delimitedList(complexSlice) | pp.delimitedList(SimpleSlice))
    + RPAREN
)
featLocation = pp.Group(SimpleSlice | complexSlice)


def parseGBLoc(s, l_, t):
    """retwingles parsed genbank location strings, assumes no joins of RC and FWD sequences"""
    strand = 1
    locationlist = []

    # see if there are any complement operators
    for entry in t[0]:
        if entry == "complement":
            strand = -1

    for entry in t[0]:
        if not isinstance(entry, str):
            locationlist.append([entry[0], entry[1]])

    # return locationlist and strand spec
    return [["location", locationlist], ["strand", strand]]


featLocation.setParseAction(parseGBLoc)

# ==== Genbank Feature Key-Value Pairs


def strip_multiline(s, l_, t):
    whitespace = re.compile("[\n]{1}[ ]+")
    return whitespace.sub(" ", t[0])


def toInt(s, l_, t):
    return int(t[0])


# Quoted KeyVal:   /key="value"
QuoteFeaturekeyval = pp.Group(
    pp.Suppress("/")
    + pp.Word(pp.alphas + pp.nums + "_-")
    + pp.Suppress("=")
    + pp.QuotedString('"', multiline=True).setParseAction(strip_multiline)
)

# UnQuoted KeyVal: /key=value  (I'm assuming it doesn't do multilines this way? wrong! ApE does store long labels this way! sigh.)
# NoQuoteFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + OneOrMore(CharsNotIn("\n")) )
keyvalspacedline = (
    pp.White(exact=21)
    + pp.CharsNotIn("/")
    + pp.OneOrMore(pp.CharsNotIn("\n"))
    + pp.LineEnd()
)
NoQuoteFeaturekeyval = pp.Group(
    pp.Suppress("/")
    + pp.Word(pp.alphas + pp.nums + "_-")
    + pp.Suppress("=")
    + pp.Combine(pp.CharsNotIn("\n") + pp.LineEnd() + pp.ZeroOrMore(keyvalspacedline))
)

# Special Case for Numerical Vals:  /bases=12  OR  /bases="12"
NumFeaturekeyval = pp.Group(
    pp.Suppress("/")
    + pp.Word(pp.alphas + pp.nums + "_-")
    + pp.Suppress("=")
    + (pp.Suppress('"') + pp.Word(pp.nums).setParseAction(toInt) + pp.Suppress('"'))
    | (pp.Word(pp.nums).setParseAction(toInt))
)

# Key Only KeyVal: /pseudo
# post-parse convert it into a pair to resemble the structure of the first three cases i.e. [pseudo, True]
FlagFeaturekeyval = pp.Group(
    pp.Suppress("/") + pp.Word(pp.alphas + pp.nums + "_-")
).setParseAction(lambda s, l_, t: [[t[0][0], True]])

Feature = pp.Group(
    pp.Word(pp.alphas + pp.nums + "_-").setParseAction(
        lambda s, l_, t: [["type", t[0]]]
    )
    + featLocation.setResultsName("location")
    + pp.OneOrMore(
        NumFeaturekeyval | QuoteFeaturekeyval | NoQuoteFeaturekeyval | FlagFeaturekeyval
    )
)

FeaturesEntry = (
    pp.Literal("FEATURES")
    + pp.Literal("Location/Qualifiers")
    + pp.Group(pp.OneOrMore(Feature)).setResultsName("features")
)

# ===============================================================================
# GenBank Sequence Parser

# sequence is just a column-spaced big table of dna nucleotides
# should it recognize full IUPAC alphabet?  NCBI uses n for unknown region
Sequence = pp.OneOrMore(
    pp.Suppress(pp.Word(pp.nums)) + pp.OneOrMore(pp.Word("ACGTacgtNn"))
)

# Group(  ) hides the setResultsName names def'd inside, such that one needs to first access this group and then access the dict of contents inside
SequenceEntry = pp.Suppress(pp.Literal("ORIGIN")) + Sequence.setParseAction(
    lambda s, l_, t: "".join(t)
).setResultsName("sequence")


# ===============================================================================
# Final GenBank Parser

# GB files with multiple records split by "//" sequence at beginning of line
GBEnd = pp.Literal("//")

# Begin w. LOCUS, slurp all entries, then stop at the end!
GB = LocusEntry + pp.OneOrMore(FeaturesEntry | SequenceEntry | GenericEntry) + GBEnd

# NCBI often returns sets of GB files
multipleGB = pp.OneOrMore(pp.Group(GB))

# ===============================================================================
# End Genbank Parser
# ===============================================================================


# ===============================================================================
# Main JSON Conversion Routine


def strip_indent(str):
    whitespace = re.compile("[\n]{1}(COMMENT){0,1}[ ]+")
    return whitespace.sub("\n", str)


def concat_dict(dlist):
    """more or less dict(list of string pairs) but merges
    vals with the same keys so no duplicates occur
    """
    newdict = {}
    for e in dlist:
        if e[0] in newdict.keys():
            newdict[e[0]] = newdict[e[0]] + strip_indent(e[1])
        else:
            newdict[e[0]] = strip_indent(e[1])
    return newdict


def toJSON(gbkstring):
    parsed = multipleGB.parseString(gbkstring)

    jseqlist = []

    for seq in parsed:
        # for item in seq.asList():
        #    print(item)

        # import sys;sys.exit(42)

        # Print to STDOUT some details (useful for long multi-record parses)
        # print(seq['name'], ":  length:", len(seq['sequence']) , " #features:" , len(seq['features'].asList()))

        # build JSON object

        nl = []
        if "features" in seq:
            for a in list(map(dict, seq["features"].asList())):
                dct = {}
                for key in a:
                    val = a[key]
                    # print(key, a[key])
                    dct[key] = a[key]
                    if isinstance(val, str):
                        dct[key] = a[key].strip()
                nl.append(dct)

        # import sys;sys.exit(42)

        # print(list(map(dict, hej))[2]["codon_start"])

        jseq = {
            "__format__": "jseq v0.1",
            "name": seq["name"],
            "size": seq["size"],
            "seqtype": seq["seqtype"],
            "divcode": seq["divcode"],
            "date": seq["date"],
            "topology": seq["topology"],
            "sequence": seq["sequence"],
            "features": nl,
            "annotations": concat_dict(seq["generics"]),
        }
        jseqlist.append(jseq)

    return jseqlist


def wrapstring(str_, rowstart, rowend, padfirst=True):
    """
    wraps the provided string in lines of length rowend-rowstart
    and padded on the left by rowstart.
    -> if padfirst is false the first line is not padded
    """
    rowlen = rowend - rowstart
    leftpad = rowstart
    wrappedstr = ""

    # no wrapping needed, single line
    if len(str_) / rowlen < 1:
        if padfirst:
            return leftpad * " " + str_ + "\n"
        else:
            return str_ + "\n"

    # multiple lines so wrap:
    for linenum in range(1 + int(len(str_) / rowlen)):
        if linenum == 0 and not padfirst:
            wrappedstr += str_[linenum * rowlen : (linenum + 1) * rowlen] + "\n"
        else:
            wrappedstr += (
                " " * leftpad + str_[linenum * rowlen : (linenum + 1) * rowlen] + "\n"
            )
    #    if str_.startswith("/translation="):
    #        print(str_)
    #        print(wrappedstr)
    #        print(".................................")
    return wrappedstr


def locstr(locs, strand):
    "genbank formatted location string, assumes no join'd combo of rev and fwd seqs"
    # slice format is like: 1..10,20..30,101..200
    locstr = ",".join(map((lambda x: str(x[0]) + ".." + str(x[1])), locs))
    if len(locs) > 1:
        locstr = "join(" + locstr + ")"
    if int(strand) == -1:
        locstr = "complement(" + locstr + ")"
    return locstr


def originstr(sequence):
    "formats dna sequence as broken, numbered lines ala genbank"
    wordlen = 10
    cols = 6
    rowlen = wordlen * cols
    outstr = ""
    for linenum in range(int(len(sequence) / rowlen) + 1):
        pos = linenum * rowlen
        # position of string for this row, then six blocks of dna
        outstr += (
            (" " * 9 + str(pos + 1))[-9:]
            + " "
            + sequence[pos : pos + 10]
            + " "
            + sequence[pos + 10 : pos + 20]
            + " "
            + sequence[pos + 20 : pos + 30]
            + " "
            + sequence[pos + 30 : pos + 40]
            + " "
            + sequence[pos + 40 : pos + 50]
            + " "
            + sequence[pos + 50 : pos + 60]
            + "\n"
        )
    return outstr


def toGB(jseq):
    "parses json jseq data and prints out ApE compatible genbank"

    # construct the LOCUS header string
    #  LOCUS format:
    #    Positions  Contents
    #    ---------  --------
    #    00:06      LOCUS
    #    06:12      spaces
    #    12:??      Locus name
    #    ??:??      space
    #    ??:40      Length of sequence, right-justified
    #    40:44      space, bp, space
    #    44:47      Blank, ss-, ds-, ms-
    #    47:54      Blank, DNA, RNA, tRNA, mRNA, uRNA, snRNA, cDNA
    #    54:55      space
    #    55:63      Blank (implies linear), linear or circular
    #    63:64      space
    #    64:67      The division code (e.g. BCT, VRL, INV)
    #    67:68      space
    #    68:79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)

    name = jseq["name"] or "default"
    size = jseq["size"] or "100"
    seqtype = jseq["seqtype"] or "DNA"
    prefix = ""
    for p in ["ds-", "ss-", "ms-"]:
        a, *b = seqtype.split(p)
        if b:
            prefix = p
            seqtype = b.pop()
            break
    prefix = prefix or "ds-"
    topology = jseq["topology"] or "linear"
    divcode = jseq["divcode"] or "   "
    date = jseq["date"] or "19-MAR-1970"

    locusstr = "LOCUS       {name:<24} {size:>4} bp {prefix}{seqtype:<4}    {topology:<8} {divcode} {date}\n".format(
        name=name,
        size=size,
        prefix=prefix,
        seqtype=seqtype,
        topology=topology,
        divcode=divcode,
        date=date,
    )

    # All these fields are left empty
    gbprops = (
        "DEFINITION  .\n"
        + "ACCESSION   \n"
        + "VERSION     \n"
        + "SOURCE      .\n"
        + "ORGANISM  .\n"
        + "COMMENT     \n"
        + "COMMENT     ApEinfo:methylated:1\n"
        + "FEATURES             Location/Qualifiers\n"
    )

    # build the feature table
    featuresstr = ""
    if "features" in jseq:
        for feat in jseq["features"]:
            fstr = (
                " " * 5
                + feat["type"]
                + " " * (16 - len(feat["type"]))
                + wrapstring(locstr(feat["location"], feat["strand"]), 21, 80, False)
            )
            for k in feat.keys():
                if k not in ["type", "location", "strand"]:
                    # ApE idiosyncrasy: don't wrap val in quotation marks
                    if k in [
                        "ApEinfo_label",
                        "ApEinfo_fwdcolor",
                        "ApEinfo_revcolor",
                        "label",
                    ]:
                        fstr += wrapstring("/" + str(k) + "=" + str(feat[k]), 21, 80)
                    # standard: wrap val in quotes
                    else:
                        fstr += wrapstring(
                            "/" + str(k) + "=" + '"' + str(feat[k]) + '"', 21, 80
                        )
            featuresstr += fstr

    # the spaced, numbered sequence
    gborigin = "ORIGIN\n" + originstr(jseq["sequence"]) + "//\n"

    return locusstr + gbprops + featuresstr + gborigin


def gbtext_clean(gbtext):
    """This function takes a string containing **one** genbank sequence
    in Genbank format and returns a named tuple containing two fields,
    the gbtext containing a string with the corrected genbank sequence and
    jseq which contains the JSON intermediate.

    Examples
    --------

    >>> s = '''LOCUS       New_DNA      3 bp    DNA   CIRCULAR SYN        19-JUN-2013
    ... DEFINITION  .
    ... ACCESSION
    ... VERSION
    ... SOURCE      .
    ...   ORGANISM  .
    ... COMMENT
    ... COMMENT     ApEinfo:methylated:1
    ... ORIGIN
    ...         1 aaa
    ... //'''
    >>> from pydna.readers import read
    >>> read(s)  # doctest: +SKIP
    ... /site-packages/Bio/GenBank/Scanner.py:1388: BiopythonParserWarning: Malformed LOCUS line found - is this correct?
    :'LOCUS       New_DNA      3 bp    DNA   CIRCULAR SYN        19-JUN-2013\\n'
      "correct?\\n:%r" % line, BiopythonParserWarning)
    Traceback (most recent call last):
      File "... /pydna/readers.py", line 48, in read
        results = results.pop()
    IndexError: pop from empty list
    <BLANKLINE>
    During handling of the above exception, another exception occurred:
    <BLANKLINE>
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "... /pydna/readers.py", line 50, in read
        raise ValueError("No sequences found in data:\\n({})".format(data[:79]))
    ValueError: No sequences found in data:
    (LOCUS       New_DNA      3 bp    DNA   CIRCULAR SYN        19-JUN-2013
    DEFINITI)
    >>> from pydna.genbankfixer import gbtext_clean
    >>> s2, j2 = gbtext_clean(s)
    >>> print(s2)
    LOCUS       New_DNA                    3 bp ds-DNA     circular SYN 19-JUN-2013
    DEFINITION  .
    ACCESSION
    VERSION
    SOURCE      .
    ORGANISM  .
    COMMENT
    COMMENT     ApEinfo:methylated:1
    FEATURES             Location/Qualifiers
    ORIGIN
            1 aaa
    //
    >>> s3 = read(s2)
    >>> s3
    Dseqrecord(o3)
    >>> print(s3.format())
    LOCUS       New_DNA                    3 bp    DNA     circular SYN 19-JUN-2013
    DEFINITION  .
    ACCESSION   New_DNA
    VERSION     New_DNA
    KEYWORDS    .
    SOURCE
      ORGANISM  .
                .
    COMMENT
                ApEinfo:methylated:1
    FEATURES             Location/Qualifiers
    ORIGIN
            1 aaa
    //"""

    jseqlist = toJSON(gbtext)
    jseq = jseqlist.pop()
    from collections import namedtuple
    from pydna._pretty import pretty_str as ps

    Result = namedtuple("Result", "gbtext jseq")
    result = Result(ps(toGB(jseq).strip()), jseq)
    return result
