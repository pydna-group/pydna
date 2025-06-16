# -*- coding: utf-8 -*-
"""
GenBank FlatFile Writer

Given a JSON representation, writes out a GenBank Flatfile (GFF) parseable by ApE

Copyright (c) 2011 Anselm Levskaya (http://anselmlevskaya.com)
Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php) license.
"""

import json
from optparse import OptionParser

# ===============================================================================
#
# This writer is a hack, much like genbank itself.
#
# Currently, it has only been verified to work with ApE and doesn't even bother
# to save a complete record of the seq-wide properties (such as comments,
# accession, etc.) as this is mainly meant for a way to inspect algorithmically
# constructed sequence from within ApE and to share sequences with others.
#
# Todo:
# write and properties into proper genbank spots, i.e. COMMENTS entries
#
# ===============================================================================


# ===============================================================================
# GenBank Formatting Helpers
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
    if len(str_) / rowlen == 0:
        if padfirst:
            return leftpad * " " + str_ + "\n"
        else:
            return str_ + "\n"

    # multiple lines so wrap:
    for linenum in range(1 + len(str_) / rowlen):
        if linenum == 0 and not padfirst:
            wrappedstr += str[linenum * rowlen : (linenum + 1) * rowlen] + "\n"
        else:
            wrappedstr += " " * leftpad + str_[linenum * rowlen : (linenum + 1) * rowlen] + "\n"
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
    for linenum in range(len(sequence) / rowlen + 1):
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


# ===============================================================================
# Main Routine


def toGB(jseqs):
    "parses json jseq data and prints out ApE compatible genbank"

    # take first jseq from parsed list
    if type(jseqs) == type([]):
        jseq = jseqs[0]
    else:
        jseq = jseqs

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
    locusstr = (
        "LOCUS"
        + " " * 6
        + jseq["name"]
        + " " * ((42 - 16) - len(jseq["name"]) - len(str(len(jseq["sequence"]))))
        + str(len(jseq["sequence"]))
        + " bp "
        + (jseq["type"] + " ")[:7]
        + " "
        + (jseq["topology"] + " ")[:8]
        + " "
        + "    "
        + jseq["date"]
        + "\n"
    )

    # fuck this noise for now
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
                    fstr += wrapstring("/" + str(k) + "=" + '"' + str(feat[k]) + '"', 21, 80)
        featuresstr += fstr

    # the spaced, numbered sequence
    gborigin = "ORIGIN\n" + originstr(jseq["sequence"]) + "//\n"

    return locusstr + gbprops + featuresstr + gborigin


# ===============================================================================
# command-line json to GB conversion
if __name__ == "__main__":
    # parse command line string
    usage = "usage: %prog [options] jseqfile_in gbfile_out"
    parser = OptionParser(usage)
    # parser.add_option("-o", "--output", dest="output",
    #                  help="output json file")
    (options, args) = parser.parse_args()

    # Load the GBK file and Build a JSON object out of the parse tree
    if len(args) > 1:
        # read in json, parse it to data
        infile = open(args[0], "r").read()
        jseq = json.loads(infile)

        # print jseq data into genbank format
        gbstr = toGB(jseq)

        # Output to new GB file
        outfile = open(args[1], "w")
        outfile.write(gbstr)
