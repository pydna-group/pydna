#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
dscode - The nucleic acid alphabet used in pydna

This file serves to define the DNA alphabet used in pydna. Each symbol usually
represents a basepair (two opposing bases in the two antiparalell DNA strands).

The alphabet is defined in this docstring which serve as the single source of
thruth for the alphabet.

A series of dictionaries:

- basepair_dict
- annealing_dict
- annealing_dict_w_holes
- complement_dict_for_dscode
- watson_tail_letter_dict
- crick_tail_letter_dict

The following bytestring translation tables are constructed from the
content in this docstring:

- complement_table_for_dscode
- dscode_to_watson_table
- dscode_to_crick_table
- dscode_to_crick_tail_table
- dscode_to_watson_tail_table
- dscode_to_full_sequence_table

The codestrings dictionary has the following keys (strings) in the order
indicated:

1. un_ambiguous_ds_dna
2. ds_rna
3. ambiguous_ds_dna
4. single_stranded_dna_rna
5. loops_dna_rna
6. mismatched_dna_rna
7. gap

Each value of the codestrings dictionary is also a string. This string has five
lines following this form:

W             1   Watson symbol
|             2   Pipe
C             3   Crick symbol
<empty line>  4
S             5   dscode symbol

W (line 1) and C (line 2) are complementary bases in a double stranded DNA
molecule and S (line 5) are the symbols of the alphabet (dscode) used to
describe the base pair above the symbol.

Line 2 must contain only the pipe character, indicating basepairing and
line 4 must be empty. The lines must be of equal length and a series ot
tests are performed to ensure the integrity of the alphabet. The string
definition follows this line and is contained in the last 13 lines of the
docstring:

un_ambiguous_ds_dna
|    ds_rna
|    |  ambiguous_ds_dna
|    |  |           single_stranded_dna_rna
|    |  |           |          loops_dna_rna
|    |  |           |          |          mismatched_dna_rna
|    |  |           |          |          |                  gap
|    |  |           |          |          |                  |
GATC UA RYMKSWHBVDN GATC••••U• -----AGCTU AAACCCGGGTTTUUUGCT •
|||| || ||||||||||| |||||||||| |||||||||| |||||||||||||||||| |
CTAG AU YRKMSWDVBHN ••••CTAG•U AGCTU----- ACGACTAGTCGTGCTUUU •

GATC UO RYMKSWHBVDN PEXIQFZJ$% 0123456789 !#{}&*()<>@:?[]=_; •

"""

from collections import namedtuple
import re as _re

# An alias for whitespace
emptyspace = chr(32)

lines = __doc__.rstrip().splitlines()[-13:]

assert not lines[-2]
assert set(lines[-4]) == {" ", "|"}

uppers = lines[-5]
pipes = lines[-4]
lowers = lines[-3]
codes = lines[-1]

assert (
    len(uppers.split())
    == len(lowers.split())
    == len(pipes.split())
    == len(codes.split())
)

names = [x.strip("| ") for x in lines[: len(codes.split())]]

codestrings = {}
for upper, pipe, lower, code, name in zip(
    uppers.split(), pipes.split(), lowers.split(), codes.split(), names
):
    codestrings[name.strip()] = f"{upper}\n{pipe}\n{lower}\n\n{code}\n".replace(
        "•", emptyspace
    )

# This string contains ascii letters not used in the alphabet
letters_not_in_dscode = "lL\"',-./\\^`|+~"

for name, codestring in codestrings.items():

    # This loops all codestrings and checks for consistency of format.
    lines = codestring.splitlines()

    assert len(lines) == 5, f'codestring["{name}"] does not have 5 lines'

    # We want the Watson, Crick and Symbol lines only
    # Second line has to be pipes ("|") and fourth has to be empty

    watsn, pipes, crick, empty, symbl = lines

    assert all(
        ln.isascii() for ln in (watsn, crick, symbl)
    ), f'codestring["{name}"] has non-ascii letters'

    assert all(
        ln.isupper() for ln in (watsn, crick, symbl) if ln.isalpha()
    ), f'codestring["{name}"] has non-uppercase letters'

    # check so that pipes contain only "|"
    assert set(pipes) == set(
        "|"
    ), f'codestring["{name}"] has non-pipe character(s) in line 2'

    # check so strings are the same length
    assert all(
        len(ln) == len(watsn) for ln in (watsn, pipes, crick, symbl)
    ), f'codestring["{name}"] has lines of unequal length'

    # These characters are not used.
    assert not any(
        [letter in letters_not_in_dscode for letter in symbl]
    ), f'codestring["{name}"] has chars outside alphabet'

"""
The `codes` dictionary is a dict of dicts containing the information of the
code strings in the form if a dict with string names as keys, each containing a
dict with this structure: (Watson letter, Crick letter): dscode symbol
"""

codes = dict()

for name, codestring in codestrings.items():

    lines = codestring.splitlines()

    watsons, _, cricks, _, symbols = lines

    # d is an alias of codes[name] used in this loop for code clarity.
    codes[name] = d = dict()

    for watson, crick, symbol in zip(watsons, cricks, symbols):
        d[watson, crick] = symbol

del d

basepair_dict = (
    codes["un_ambiguous_ds_dna"]
    | codes["ambiguous_ds_dna"]
    | codes["ds_rna"]
    | codes["single_stranded_dna_rna"]
    # | codes["mismatched_dna_rna"]
    # | codes["loops_dna_rna"]
    | codes["gap"]
)

annealing_dict = dict()

"""
The annealing_dict_of_str is constructed below. It contains the information needed
to tell if two DNA fragments (like a and b below) can anneal. This of cource only concerns
single stranded regions.

The dict has the form (x, y): s
Where x and y are bases in a and b and the symbol s is the resulting symbol for the base pair
that is formed. The letters x and y are from the values of the codes["single_stranded_dna_rna"]
dictionary.

For, example: One key-value pair is ('P', 'Q'): 'G' which matches the first
of the four new base pairings formed between a and b in the example below.


(a)
gggPEXI    (dscode for a)

gggGATC
ccc
       aaa (b)
   CTAGttt

   QFZJaaa (dscode for b)


gggGATCaaa (annealing product between a and b)
cccCTAGttt

This loops through the base pairs where the upper or lower
positions are empty. (w, c), s would be ("G", " "), "P"
in the first iteration.
"""

temp = codes["un_ambiguous_ds_dna"] | codes["ds_rna"]

# Alias to make the code below more readable.
d = codes["single_stranded_dna_rna"]

for (x, y), symbol in d.items():
    if y == emptyspace:
        other = next(b for a, b in temp if a == x)
        symbol_other = d[emptyspace, other]
        annealing_dict[symbol, symbol_other] = temp[x, other]
        annealing_dict[symbol_other, symbol] = temp[x, other]
    elif x == emptyspace:
        other = next(a for a, b in temp if b == y)
        symbol_other = d[other, emptyspace]
        annealing_dict[symbol, symbol_other] = temp[other, y]
        annealing_dict[symbol_other, symbol] = temp[other, y]
    else:
        raise ValueError("This should not happen")

del d, temp

temp = {}

for (x, y), symbol in annealing_dict.items():

    temp[x, emptyspace] = x
    temp[emptyspace, y] = y

annealing_dict_w_holes = annealing_dict | temp

del temp

"""
A collection of translation tables are a practical way to obtain Watson and Crick
from dscode or the reverse complement strands when needed.

These are meant to be used by the bytes.translate method.


The translation table "complement_table_for_dscode" is used to obtain the
complement of a DNA sequence in dscode format.
"""

complement_dict_for_dscode = {
    s: basepair_dict[c, w] for (w, c), s in basepair_dict.items()
}

from_letters = "".join(complement_dict_for_dscode.keys())
to_letters = "".join(complement_dict_for_dscode.values())

from_letters += from_letters.lower()
to_letters += to_letters.lower()

complement_table_for_dscode = bytes.maketrans(
    from_letters.encode("ascii"), to_letters.encode("ascii")
)

"""
dscode_to_watson_table and dscode_to_crick_table are used to obtain the Watson
and (reverse) Crick strands from dscode. Four extra letters are added to the
table and used in the pydna.dseq.representation function. These are used
to add range indicators ("..") in the watson or crick strings for
representation of long sequences.

The four letters are placeholder1, placeholder2, interval, empty_bs
"""

dscode_sense = ""
dscode_compl = ""
watson = ""
crick = ""
dscode_sense_lower = ""
dscode_compl_lower = ""
watson_lower = ""
crick_lower = ""

for (w, c), dscode in basepair_dict.items():
    dscode_sense += dscode
    dscode_compl += basepair_dict[c, w]
    watson += w
    crick += c
    dscode_lower = dscode.lower()
    if dscode_lower in dscode_sense:
        continue
    dscode_sense_lower += dscode_lower
    watson_lower += w.lower()
    crick_lower += c.lower()
    dscode_compl_lower += dscode_compl.lower()

# dscode_sense += dscode_sense.lower()
# dscode_compl += dscode_compl.lower()
# watson += watson.lower()
# crick += crick.lower()

placeholder1 = "~"
placeholder2 = "+"
interval = "."

assert placeholder1 in letters_not_in_dscode
assert placeholder2 in letters_not_in_dscode
assert interval in letters_not_in_dscode

dscode_to_watson_table = bytes.maketrans(
    (dscode_sense + dscode_sense_lower + placeholder1 + placeholder2).encode("ascii"),
    (watson + watson_lower + emptyspace + interval).encode("ascii"),
)

dscode_to_crick_table = bytes.maketrans(
    (dscode_sense + dscode_sense_lower + placeholder1 + placeholder2).encode("ascii"),
    (crick + crick_lower + interval + emptyspace).encode("ascii"),
)


watson_tail_letter_dict = {
    w: s for (w, c), s in codes["single_stranded_dna_rna"].items() if c.isspace()
}

from_letters = "".join(watson_tail_letter_dict.keys())
to_letters = "".join(watson_tail_letter_dict.values())

from_letters += from_letters.lower()
to_letters += to_letters.lower()

dscode_to_crick_tail_table = bytes.maketrans(
    from_letters.encode("ascii"), to_letters.encode("ascii")
)

from_letters_full = five_prime_ss_letters = to_letters
to_letters_full = from_letters


d = codes["single_stranded_dna_rna"]

crick_tail_letter_dict = {
    complement_dict_for_dscode[c]: s for (w, c), s in d.items() if w.isspace()
}

del d

from_letters = "".join(crick_tail_letter_dict.keys())
to_letters = "".join(crick_tail_letter_dict.values())

from_letters += from_letters.lower()
to_letters += to_letters.lower()

dscode_to_watson_tail_table = bytes.maketrans(
    from_letters.encode("ascii"), to_letters.encode("ascii")
)

three_prime_ss_letters = to_letters
from_letters_full += to_letters
to_letters_full += from_letters

dscode_to_full_sequence_table = bytes.maketrans(
    from_letters_full.encode("ascii"), to_letters_full.encode("ascii")
)


# This loop adds upper and lower case symbols
mixed_case_dict = {}

for (x, y), symbol in basepair_dict.items():
    mixed_case_dict[x.lower(), y.lower()] = symbol.lower()
    mixed_case_dict[x.lower(), y.upper()] = symbol.lower()
    mixed_case_dict[x.upper(), y.lower()] = symbol.upper()

    if x == emptyspace:
        mixed_case_dict[x, y.lower()] = symbol.lower()
        mixed_case_dict[x, y.upper()] = symbol.upper()
    if y == emptyspace:
        mixed_case_dict[x.lower(), y] = symbol.lower()
        mixed_case_dict[x.upper(), y] = symbol.upper()

# Add mixed case entries to the dict
basepair_dict.update(mixed_case_dict)

# This loop adds upper and lower case symbols
mixed_case_dict = {}

for (x, y), symbol in annealing_dict.items():
    mixed_case_dict[x.lower(), y.lower()] = symbol.lower()
    mixed_case_dict[x.lower(), y.upper()] = symbol.lower()
    mixed_case_dict[x.upper(), y.lower()] = symbol.upper()
# Add mixed case entries to the dict
annealing_dict.update(mixed_case_dict)

ds_letters = (
    "".join(codes["un_ambiguous_ds_dna"].values())
    + "".join(codes["ds_rna"].values())
    + "".join(codes["ambiguous_ds_dna"].values())
)

ss_letters_watson = "".join(
    s for (w, c), s in codes["single_stranded_dna_rna"].items() if c == emptyspace
)
ss_letters_crick = "".join(
    s for (w, c), s in codes["single_stranded_dna_rna"].items() if w == emptyspace
)

ds_letters += ds_letters.lower()
ss_letters_watson += ss_letters_watson.lower()
ss_letters_crick += ss_letters_crick.lower()


"""
The dict of regexes below cover IUPAC Ambiguity Code complements
and is used in the amplify module.
"""
iupac_compl_regex = {
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


# This loop adds upper and lower case symbols
# mixed_case_dict = {}

for (x, y), symbol in annealing_dict_w_holes.items():
    mixed_case_dict[x.lower(), y.lower()] = symbol.lower()
    mixed_case_dict[x.lower(), y.upper()] = symbol.lower()
    mixed_case_dict[x.upper(), y.lower()] = symbol.upper()
# Add mixed case entries to the dict
annealing_dict_w_holes.update(mixed_case_dict)


def get_parts(datastring: str) -> namedtuple:
    """
    A namedtuple containing the parts of a dsDNA sequence.

    The datastring should contain a string with dscode symbols.
    A regex is used to capture the single stranded regions at the ends as
    well as the regiond in the middle.

    The figure below numbers the regex capture groups and what they capture
    as well as the namedtuple field name.

    ::

         group 0 "sticky_left5"
         |
         |      group 3"sticky_right5"
         |      |
        ---    ---
        GGGATCC
           TAGGTCA
           ----
             |
             group 2 "middle"



         group 1 "sticky_left3"
         |
         |      group 4 "sticky_right3"
         |      |
        ---    ---
           ATCCAGT
        CCCTAGG
           ----
             |
             group 2 "middle"



           group 5 "single_watson" (only an upper strand)
           |
        -------
        ATCCAGT
        |||||||



           group 6 "single_crick" (only a lower strand)
           |
        -------

        |||||||
        CCCTAGG


    Up to seven groups (0..6) are captured, but some are mutually exclusive
    which means that one of them is an empty string:

    0 or 1, not both, a DNA fragment has either 5' or 3' sticky end.

    2 or 5 or 6, a DNA molecule has a ds region or is single stranded.

    3 or 4, not both, either 5' or 3' sticky end.

    Note that internal single stranded regions are not identified and will
    be contained in the middle part if they are present.

    Parameters
    ----------
    datastring : str
        A string with dscode.

    Returns
    -------
    namedtuple
        Seven string fields describing the DNA molecule.
        fragment(sticky_left5='', sticky_left3='',
                 middle='',
                 sticky_right3='', sticky_right5='',
                 single_watson='', single_crick='')

    """

    m = _re.match(
        f"([{ss_letters_watson}]*)"  # capture group 0 ssDNA in watson strand
        f"([{ss_letters_crick}]*)"  # "             1 ssDNA in crick strand
        f"(?=[{ds_letters}])"  # positive lookahead for dsDNA, no capture
        "(.*)"  # capture group 2 everything in the middle
        f"(?<=[{ds_letters}])"  # positive look behind for dsDNA, no capture
        f"([{ss_letters_watson}]*)"  # capture group 3 ssDNA in watson strand
        f"([{ss_letters_crick}]*)|"  # "             4 ssDNA in crick strand
        f"([{ss_letters_watson}]+)|"  # "             5 if data contains only upper strand
        f"([{ss_letters_crick}]+)",  # "             6 if data contains only lower strand
        datastring,
    )

    result = m.groups() if m else (None, None, None, None, None, None, None)

    result = ["" if e is None else e for e in result]

    field_names = (
        "sticky_left5",
        "sticky_left3",
        "middle",
        "sticky_right3",
        "sticky_right5",
        "single_watson",
        "single_crick",
    )

    fragment = namedtuple("fragment", field_names)

    return fragment(*result)


def dsbreaks(data: str):

    wl = _re.escape(five_prime_ss_letters)
    cl = _re.escape(three_prime_ss_letters)

    breaks = []
    regex = (
        "(.{0,3})"  # context if present
        f"([{wl}][{cl}]|[{cl}][{wl}])"  # ss chars next to each other
        "(.{0,3})"  # context if present
    )
    for mobj in _re.finditer(regex, data):
        chunk = mobj.group()
        w, c = representation_tuple(chunk)
        breaks.append(f"[{mobj.start()}:{mobj.end()}]\n{w}\n{c}\n")
    return breaks


def representation_tuple(
    datastring: str = "", length_limit_for_repr: int = 30, chunk: int = 4
):
    """
    Two line string representation of a sequence of dscode symbols.

    See pydna.alphabet module for the definition of the pydna dscode
    alphabet. The dscode has a symbol (ascii) character for base pairs
    and single stranded DNA.

    This function is used by the Dseq.__repr__() method.

    Parameters
    ----------
    data : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    str
        A two line string containing The Watson and Crick strands.

    """

    (
        sticky_left5,
        sticky_left3,
        middle,
        sticky_right5,
        sticky_right3,
        single_watson,
        single_crick,
    ) = get_parts(datastring)

    if len(datastring) > length_limit_for_repr:
        """
        We need to shorten the repr if the sequence is longer than
        limit imposed by length_limit_for_repr.

        The representation has three parts, so we divide by three for each part.

        Long DNA strands are interrupted by interval notation, like agc..att
        where the two dots indicate intervening hidden sequence.


        Dseq(-71)
        GAAA..AATCaaaa..aaaa
                  tttt..ttttCTAA..AAAG

        placeholder1, placeholder2 are two letters that are replaced by
        interval characters in the upper or lower strands by the translation
        """

        part_limit = length_limit_for_repr // 3

        if len(sticky_left5) > part_limit:
            sticky_left5 = (
                sticky_left5[:chunk] + placeholder2 * 2 + sticky_left5[-chunk:]
            )

        if len(sticky_left3) > part_limit:
            sticky_left3 = (
                sticky_left3[:chunk] + placeholder1 * 2 + sticky_left3[-chunk:]
            )

        if len(middle) > part_limit:
            middle = middle[:4] + interval * 2 + middle[-4:]

        if len(sticky_right5) > part_limit:
            sticky_right5 = (
                sticky_right5[:chunk] + placeholder2 * 2 + sticky_right5[-chunk:]
            )

        if len(sticky_right3) > part_limit:
            sticky_right3 = (
                sticky_right3[:chunk] + placeholder1 * 2 + sticky_right3[-chunk:]
            )

    """
    The processed string contains
    """
    processed_dscode = (sticky_left5 or sticky_left3) + middle + (
        sticky_right5 or sticky_right3
    ) or single_watson + single_crick

    watson = processed_dscode.translate(dscode_to_watson_table).rstrip()
    crick = processed_dscode.translate(dscode_to_crick_table).rstrip()

    return watson, crick
