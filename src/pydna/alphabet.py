#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Six multiline strings are defined in this file.

codestrings["un_ambiguous_ds_dna"]
codestrings["ambiguous_ds_dna"]
codestrings["ds_rna"]
codestrings["single_stranded_dna_rna"]
codestrings["mismatched_dna_rna"]
codestrings["loops_dna_rna"]

Each string has five lines and describe the DNA alphabet
used in Pydna in this form:

W             1
|             2
C             3
<empty line>  4
S             5

W (line 1) and C (line 2) are complementary bases in a double stranded DNA molecule and S (line 5) are
the symbols of the alphabet used to describe the base pair above the symbol.

"""

emptyspace = chr(32)

codestrings = dict()


codestrings[
    "un_ambiguous_ds_dna"
] = """\
GATC
||||
CTAG

GATC
"""

codestrings[
    "ambiguous_ds_dna"
] = """\
RYMKSWHBVDN
|||||||||||
YRKMSWDVBHN

RYMKSWHBVDN
"""

codestrings[
    "ds_rna"
] = """\
UA
||
AU

UO
"""

codestrings["single_stranded_dna_rna"] = (
    """\
GATC....U.
||||||||||
....CTAG.U

PEXIQFZJ$%
""".replace(
        ".", emptyspace
    )
)

codestrings[
    "mismatched_dna_rna"
] = """\
AAACCCGGGTTTUUUGCT
||||||||||||||||||
ACGACTAGTCGTGCTUUU

!#{}&*()<>@:?[]=_;
"""

codestrings[
    "loops_dna_rna"
] = """\
-----AGCTU
||||||||||
AGCTU-----

0123456789
"""

keys = set(
    (
        "un_ambiguous_ds_dna",
        "ambiguous_ds_dna",
        "ds_rna",
        "single_stranded_dna_rna",
        "mismatched_dna_rna",
        "loops_dna_rna",
    )
)

assert set(codestrings.keys()) == keys

not_dscode = "lL\"',-./\\^`|+~"

for name, codestring in codestrings.items():

    # This loops all codestrings and checks for consistency of format.
    lines = codestring.splitlines()

    assert len(lines) == 5

    # We want the Watson, Crick and Symbol lines only
    # Second line has to be pipes ("|") and fourth has to be empty

    watsn, pipes, crick, empty, symbl = lines

    assert all(ln.isascii() for ln in (watsn, crick, symbl))

    assert all(ln.isupper() for ln in (watsn, crick, symbl) if ln.isalpha())

    # check so that pipes contain only "|"
    assert set(pipes) == set("|")

    # check so strings are the same length
    assert all(len(ln) == len(watsn) for ln in (watsn, pipes, crick, symbl))

    # These characters are not used.
    assert not any([letter in not_dscode for letter in symbl])


codes = dict()

for name, codestring in codestrings.items():

    lines = codestring.splitlines()

    watsons, _, cricks, _, symbols = lines

    codes[name] = dct = dict()

    for watson, crick, symbol in zip(watsons, cricks, symbols):
        if watson == emptyspace:
            dct[watson, crick.lower()] = symbol.lower()
            dct[watson, crick.upper()] = symbol.upper()
        else:
            dct[watson.upper(), crick.upper()] = symbol.upper()
            dct[watson.upper(), crick.lower()] = symbol.upper()
            dct[watson.lower(), crick.upper()] = symbol.lower()
            dct[watson.lower(), crick.lower()] = symbol.lower()


bp_dict_str = (
    codes["un_ambiguous_ds_dna"]
    | codes["ambiguous_ds_dna"]
    | codes["ds_rna"]
    | codes["single_stranded_dna_rna"]
    # | codes["mismatched_dna_rna"]
    # | codes["loops_dna_rna"]
)

bp_dict = {
    (w.encode("ascii"), c.encode("ascii")): s.encode("ascii")
    for (w, c), s in bp_dict_str.items()
}

temp = codes["un_ambiguous_ds_dna"] | codes["ds_rna"]


annealing_dict_str = dict()

# The annealing_dict_str is constructed below. This dict contains the information needed
# to tell if two DNA fragments (like a and b below) can anneal. The dict has the form (x, y): s
# Where x and y are bases in a and b and the symbol s is the resulting symbol for the base pair
# that is formed. One element in the dict is ('P', 'Q'): 'G' which matches the first
# of the four new base pairings formed between a and b in the example below.
#
#
# (a)
# gggPEXI    (dscode for a)
#
# gggGATC
# ccc
#        aaa (b)
#    CTAGttt
#
#    QFZJaaa (dscode for b)
#
#
# gggGATCaaa (annealing product between a and b)
# cccCTAGttt
#
# This loops through the base pairs where the upper or lower
# positions are empty. (w, c), s would be ("G", " "), "P"
# in the first iteration.
#

d = codes["single_stranded_dna_rna"]  # Alias to make the code below more readable.

for (x, y), symbol in d.items():
    if y == emptyspace:
        other = next(b for a, b in temp if a == x)
        symbol_other = d[emptyspace, other]
        annealing_dict_str[symbol, symbol_other] = temp[x, other]
        annealing_dict_str[symbol_other, symbol] = temp[x, other]
    elif x == emptyspace:
        other = next(a for a, b in temp if b == y)
        symbol_other = d[other, emptyspace]
        annealing_dict_str[symbol, symbol_other] = temp[other, y]
        annealing_dict_str[symbol_other, symbol] = temp[other, y]
    else:
        raise ValueError("This should not happen")

del d

mixed_case_dict = (
    dict()
)  # This dict will contain upper and lower case symbols annealing_dict_str

for (x, y), symbol in annealing_dict_str.items():
    mixed_case_dict[x.upper(), y.lower()] = symbol.upper()
    mixed_case_dict[x.lower(), y.upper()] = symbol.lower()
    mixed_case_dict[x.lower(), y.lower()] = symbol.lower()

annealing_dict_str = (
    annealing_dict_str | mixed_case_dict
)  # Add mixed case entries to the dict

# A bytestr version of the annealing_dict_str
annealing_dict = {
    (x.encode("ascii"), y.encode("ascii")): s.encode("ascii")
    for (x, y), s in annealing_dict_str.items()
}

dscode_sense = []
dscode_compl = []
watson = []
crick = []

for (w, c), s in bp_dict.items():

    if w.isupper() and c.islower() or w.islower() and c.isupper():
        continue

    dscode_sense.append(s)
    dscode_compl.append(bp_dict[c, w])
    watson.append(w)
    crick.append(c)

complement_table_dscode = bytes.maketrans(
    b"".join(dscode_sense), b"".join(dscode_compl)
)

placeholder1, placeholder2, interval, empty_bs = (
    b"~",
    b"+",
    b".",
    emptyspace.encode("ascii"),
)

for bstring in placeholder1, placeholder2, interval:
    assert all(letter in not_dscode.encode("ascii") for letter in bstring)

dscode_to_watson_table = bytes.maketrans(
    b"".join(dscode_sense) + placeholder1 + placeholder2,
    b"".join(watson) + empty_bs + interval,
)

dscode_to_crick_table = bytes.maketrans(
    b"".join(dscode_sense) + placeholder1 + placeholder2,
    b"".join(crick) + interval + empty_bs,
)


watson_tail_letter_dict = {
    (w.encode("ascii")): s.encode("ascii")
    for (w, c), s in codes["single_stranded_dna_rna"].items()
    if c.isspace()
}

from_letters = b"".join(watson_tail_letter_dict.keys())

to_letters = b"".join(watson_tail_letter_dict.values())

dscode_to_crick_tail_table = bytes.maketrans(from_letters, to_letters)
# dscode_to_crick_tail_table = bytes.maketrans(b"GATCgatc", b"PEXIpexi")


crick_tail_letter_dict = {
    (c.encode("ascii")): s.encode("ascii")
    for (w, c), s in codes["single_stranded_dna_rna"].items()
    if w.isspace()
}

from_letters = b"".join(crick_tail_letter_dict.keys())

to_letters = b"".join(crick_tail_letter_dict.values())

dscode_to_watson_tail_table = bytes.maketrans(from_letters, to_letters)
dscode_to_watson_tail_table = bytes.maketrans(b"GATCgatc", b"QFZJqfzj")


dscode_to_to_full_sequence_table = bytes.maketrans(
    b"PEXIpexiQFZJqfzj", b"GATCgatcGATCgatc"
)


iupac_compl_regex = {  # IUPAC Ambiguity Code complements
    "A": "(?:T|U)",
    "C": "(?:G)",
    "G": "(?:C)",
    "T": "(?:A)",
    "U": "(?:A)",
    "R": "(?:T|C|Y)",
    "Y": "(?:G|A|R)",
    "S": "(?:G|C|S)",
    "W": "(?:A|T|W)",
    "K": "(?:C|AM)",
    "M": "(?:T|G|K)",
    "B": "(?:C|G|A|V)",
    "D": "(?:A|C|T|H)",
    "H": "(?:A|G|T|D)",
    "V": "(?:T|C|G|B)",
    "N": "(?:A|G|C|T|N)",
}
