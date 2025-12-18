#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
dscode - The nucleic acid alphabet used in pydna

This file serves to define dscode, the DNA alphabet used in pydna.
Each symbol represents a basepair (two opposing bases in the two antiparalell
DNA strands).

The alphabet is defined in the end of this docstring which serve as the single
source of thruth. The alphabet is used to construct the codestrings dictionary
with has the following keys (strings) in the order indicated:

1. un_ambiguous_ds_dna
2. ds_rna
3. ambiguous_ds_dna
4. single_stranded_dna_rna
5. loops_dna_rna
6. mismatched_dna_rna
7. gap

Each value of the codestrings dictionary is a multiline string. This string
has five lines following this form:

::

    W             1   Watson symbol
    |             2   Pipe
    C             3   Crick symbol
    <empty line>  4
    S             5   dscode symbol

W (line 1) and C (line 2) are complementary bases in a double stranded DNA
molecule and S (line 5) are the symbols of the alphabet used to
describe the base pair above the symbol.

Line 2 must contain only the pipe character, indicating basepairing and
line 4 must be empty. The lines must be of equal length and a series ot
tests are performed to ensure the integrity of the alphabet.

The string definition as well as the keys for the codestrings dict follow this
line and is contained in the last 13 lines of the docstring:

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
import re
from dataclasses import dataclass

__all__ = [
    # Core alphabet dictionaries
    "basepair_dict",
    "annealing_dict",
    "annealing_dict_w_holes",
    "complement_dict_for_dscode",
    # Translation tables (str.translate, bytes.translate)
    "complement_table_for_dscode",
    "dscode_to_watson_table",
    "dscode_to_crick_table",
    "dscode_to_watson_tail_table",
    "dscode_to_crick_tail_table",
    "dscode_to_full_sequence_table",
    # Alphabet subsets
    "ds_letters",
    "ss_letters_watson",
    "ss_letters_crick",
    # Regex helpers and factories
    "iupac_compl_regex",
    "regex_ss_melt_factory",
    "regex_ds_melt_factory",
    # Data structures
    "DseqParts",
    # Public helper functions
    "get_parts",
    "dsbreaks",
    "representation_tuple",
    "anneal_strands",
]


# An alias for whitespace
emptyspace = chr(32)

# ============================================================================
# Alphabet definition extracted from module docstring
# ============================================================================

lines = __doc__.rstrip().splitlines()[-13:]  # last 13 docstring lines are read

assert not lines[-2]  # line 4 has to be empty
assert set(lines[-4]) == {" ", "|"}  # line 2 has to have pipes only.

uppers = lines[-5]
pipes = lines[-4]
lowers = lines[-3]
dscode = lines[-1]

# Make sure all lineas are equal in length
assert (
    len(uppers.split())
    == len(lowers.split())
    == len(pipes.split())
    == len(dscode.split())
)

# Extract the keys from the docstring
names = [x.strip("| ") for x in lines[: len(dscode.split())]]

# ============================================================================
# Construct the codestrings dict
# ============================================================================

codestrings = {}

for upper, pipe, lower, code, name in zip(
    uppers.split(), pipes.split(), lowers.split(), dscode.split(), names
):
    codestrings[name.strip()] = f"{upper}\n{pipe}\n{lower}\n\n{code}\n".replace(
        "•", emptyspace
    )


# ============================================================================
# Define ascii letters not used in the alphabet
# ============================================================================

letters_not_in_dscode = "lL\"',-./\\^`|+~"


# ============================================================================
# for loop below carries out a series of consistency checks
# ============================================================================

for name, codestring in codestrings.items():

    lines = codestring.splitlines()

    assert len(lines) == 5, f'codestring["{name}"] does not have 5 lines'

    # We want the Watson, Crick and Symbol lines only
    # Second line has to be pipes ("|") and fourth has to be empty

    watsn, pipes, crick, empty, symbl = lines

    # Check so that all letters are ascii symbols.
    assert all(
        ln.isascii() for ln in (watsn, crick, symbl)
    ), f'codestring["{name}"] has non-ascii letters'

    # Verify so that all chars that have uppercase are uppercase.
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

    # Check that the the letters in the letters_not_in_dscode string
    # are not used.
    assert not any(
        [letter in letters_not_in_dscode for letter in symbl]
    ), f'codestring["{name}"] has chars outside alphabet'


# ============================================================================
# The `codes` dictionary is a dict of dicts containing the information of the
# code strings in the form if a dict with string names as keys, each containing
# a {tuple: string} dict with this structure:
#
#                                  (Watson letter, Crick letter): dscode symbol
# ============================================================================

codes = dict()

for name, codestring in codestrings.items():

    lines = codestring.splitlines()

    watsons, _, cricks, _, symbols = lines

    # d is an alias of codes[name] used in this loop for code clarity.
    codes[name] = d = dict()

    for watson, crick, symbol in zip(watsons, cricks, symbols):
        d[watson, crick] = symbol

del d  # delete alias

# ============================================================================
# The `basepair_dict` dictionary is a merge of a subset of the `codes`dict.
# ============================================================================

basepair_dict = (
    codes["un_ambiguous_ds_dna"]
    | codes["ambiguous_ds_dna"]
    | codes["ds_rna"]
    | codes["single_stranded_dna_rna"]
    # | codes["mismatched_dna_rna"]
    # | codes["loops_dna_rna"]
    | codes["gap"]
)


# ============================================================================
# The `annealing_dict` dictionary contain letters for single stranded
# DNA and their dscode after annealing
# ============================================================================

# The annealing_dict_of_str is constructed below. It contains the information
# needed to tell if two DNA fragments (like a and b below) can anneal.

# This of course only concerns single stranded regions.

# The dict has the form (x, y): s

# Where x and y are bases in a and b and the symbol s is the resulting dscode
# symbol for the base pair that is formed.

# The letters x and y are from the values in the
# codes["single_stranded_dna_rna"] dictionary.

# For, example: One key-value pair is ('P', 'Q'): 'G' which matches the first
# of the four new base pairings formed between a and b in the example below.

#   (a)
#   gggPEXI    (dscode for a)

#   gggGATC
#   ccc
#           aaa (b)
#       CTAGttt

#       QFZJaaa (dscode for b)


#   gggGATCaaa  (annealing product between a and b)
#   cccCTAGttt

# This loops through the base pairs where the upper or lower
# positions are empty. (w, c), s would be ("G", " "), "P"
# in the first iteration.

annealing_dict = dict()

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

# ============================================================================
# The `annealing_dict_w_holes`contains the `annealing_dict`
# and additional key pairs where one position is empty
# ============================================================================

temp = {}

for (x, y), symbol in annealing_dict.items():

    temp[x, emptyspace] = x
    temp[emptyspace, y] = y

annealing_dict_w_holes = annealing_dict | temp

del temp


# ============================================================================
# translation tables
# ============================================================================

# A collection of translation tables are a practical way to obtain Watson and Crick
# from dscode or the reverse complement strands when needed.

# These are meant to be used by the str.translate or bytes.translate methods.


# ============================================================================
# The translation table "complement_table_for_dscode" is used to obtain the
# complement of a DNA sequence in dscode format.
# ============================================================================

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


# ============================================================================
# dscode_to_watson_table and dscode_to_crick_table
# ============================================================================

# dscode_to_watson_table and dscode_to_crick_table are used to obtain the Watson
# and (reverse) Crick strands from dscode.

# Three extra letters (placeholder1, placeholder2, interval) are added to the
# table and used in the representation_tuple function to
# add range indicators ("..") in the watson or crick strings for
# representation of long sequences.

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


# ============================================================================
# dscode_to_watson_tail_table
# ============================================================================


watson_tail_letter_dict = {
    w: s for (w, c), s in codes["single_stranded_dna_rna"].items() if c.isspace()
}

from_letters = "".join(watson_tail_letter_dict.keys())
to_letters = "".join(watson_tail_letter_dict.values())

from_letters += from_letters.lower()
to_letters += to_letters.lower()

dscode_to_watson_tail_table = bytes.maketrans(
    from_letters.encode("ascii"), to_letters.encode("ascii")
)

from_letters_full = five_prime_ss_letters = to_letters
to_letters_full = from_letters

# ============================================================================
# dscode_to_crick_tail_table
# ============================================================================

crick_tail_letter_dict = {
    complement_dict_for_dscode[c]: s
    for (w, c), s in codes["single_stranded_dna_rna"].items()
    if w.isspace()
}

from_letters = "".join(crick_tail_letter_dict.keys())
to_letters = "".join(crick_tail_letter_dict.values())

from_letters += from_letters.lower()
to_letters += to_letters.lower()

dscode_to_crick_tail_table = bytes.maketrans(
    from_letters.encode("ascii"), to_letters.encode("ascii")
)

three_prime_ss_letters = to_letters
from_letters_full += to_letters
to_letters_full += from_letters


# ============================================================================
# dscode_to_full_sequence_table
# ============================================================================


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

mixed_case_dict = {}

# This loop adds upper and lower case symbols
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


# ============================================================================
# iupac_compl_regex dict of regexes below cover IUPAC Ambiguity Code
# complements and is used in the amplify module.
# ============================================================================

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

# ============================================================================
# DseqParts dataclass
# ============================================================================


@dataclass
class DseqParts:
    sticky_left5: str
    sticky_left3: str
    middle: str
    sticky_right3: str
    sticky_right5: str
    single_watson: str
    single_crick: str

    def __iter__(self):
        """
        Allow unpacking DseqParts instances.
        >>> from pydna.alphabet import get_parts
        >>> sticky_left5, sticky_left3, middle, sticky_right3, sticky_right5, single_watson, single_crick = get_parts("eeATCGuggCCGgg")
        >>> sticky_left5
        'ee'
        >>> middle
        'ATCGuggCCGgg'
        """
        return iter(
            (
                self.sticky_left5,
                self.sticky_left3,
                self.middle,
                self.sticky_right3,
                self.sticky_right5,
                self.single_watson,
                self.single_crick,
            )
        )

    def __getitem__(self, index: int) -> str:
        """
        Allow indexing DseqParts instances.
        >>> from pydna.alphabet import get_parts
        >>> parts = get_parts("eeATCGuggCCGgg")
        >>> parts[0]
        'ee'
        >>> parts[2]
        'ATCGuggCCGgg'
        """
        return tuple(self)[index]


def get_parts(datastring: str) -> DseqParts:
    """
    Returns a DseqParts instance containing the parts of a dsDNA sequence.

    The datastring argument should contain a string with dscode symbols.

    A regular expression is used to capture the single stranded regions at
    the ends as well as the ds region in the middle, if any.

    The figure below numbers the regex capture groups and what they capture
    as well as the DseqParts instance field name for each group.

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

    Examples
    --------
    >>>

    Up to seven groups (0..6) are captured.s ome are mutually exclusive
    which means that one of them is an empty string:

    0 or 1, not both, a DNA fragment has either 5' or 3' sticky end.

    2 or 5 or 6, a DNA molecule has a ds region or is entirely single stranded.

    3 or 4, not both, either 5' or 3' sticky end.

    Note that internal single stranded regions are not identified and will
    be contained in the middle part if they are present.

    Parameters
    ----------
    datastring : str
        A string with dscode.

    Returns
    -------
    DseqParts
        Seven string fields describing the DNA molecule.
        DseqParts(sticky_left5='', sticky_left3='',
                  middle='',
                  sticky_right3='', sticky_right5='',
                  single_watson='', single_crick='')

    """

    m = re.match(
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

    return DseqParts(
        sticky_left5=result[0],
        sticky_left3=result[1],
        middle=result[2],
        sticky_right3=result[3],
        sticky_right5=result[4],
        single_watson=result[5],
        single_crick=result[6],
    )


def dsbreaks(datastring: str) -> list[str]:
    """
    Find double strand breaks in DNA in dscode format.

    An empty watson position next to an empty crick position in the dsDNA
    leads to a discontinuous DNA. This function is used to show breaks in
    DNA in Dseq.__init__.

    >>> from pydna.alphabet import dsbreaks
    >>> x, = dsbreaks("GATPFTAA")
    >>> print(x)
    [0:8]
    GATG TAA
    CTA TATT
    >>> dsbreaks("GATC")
    []

    Parameters
    ----------
    data : str
        A string representing DNA in dscode format.

    Returns
    -------
    list[str]
        A list of 3-line

    """

    wl = re.escape(five_prime_ss_letters)
    cl = re.escape(three_prime_ss_letters)

    breaks = []
    regex = (
        "(.{0,3})"  # return context if present.
        f"([{wl}][{cl}]|[{cl}][{wl}])"  # find adjacent single strand chars.
        "(.{0,3})"  # return context if present.
    )
    for mobj in re.finditer(regex, datastring):
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

    # The processed string that will be used to
    # obtain a watson and crick strand
    processed_dscode = (sticky_left5 or sticky_left3) + middle + (
        sticky_right5 or sticky_right3
    ) or single_watson + single_crick

    watson = processed_dscode.translate(dscode_to_watson_table).rstrip()
    crick = processed_dscode.translate(dscode_to_crick_table).rstrip()

    return watson, crick


def regex_ss_melt_factory(length: int) -> re.Pattern:
    """
    A regular expression for finding double-stranded regions flanked by single-stranded DNA
    that can be melted to shed a single-stranded fragment.

    This function returns a regular expression that finds double-stranded regions
    (of length <= length) that are flanked by single-stranded regions on the same
    side in dscode format. These regions are useful to identify as potential melt
    sites, since melting them leads to the shedding of a single-stranded fragment.

    The regular expression finds double stranded patches flanked by empty
    positions on the same side (see figure below). Melting of this kind of
    sites leads to the shedding of a single stranded fragment.

    ::

        GFTTAJA   <-- dscode representing the ds DNA below.

        G TTA A   <-- "TTA" is found by the regex for length <= 3
        CTAATGT


    Examples
    --------
    >>> from pydna.dseq import Dseq
    >>> regex = regex_ss_melt_factory(3)
    >>> s = Dseq("GFTTAJA")
    >>> s
    Dseq(-7)
    G TTA A
    CTAATGT
    >>> mobj = regex.search(s._data)
    >>> mobj.groupdict()
    {'watson': b'TTA', 'crick': None}


    Parameters
    ----------
    length : int
        Max length of double stranded region flanked by single stranded
        regions.

    Returns
    -------
    TYPE
        regular expression object.

    """

    regex = (
        f"(?P<watson>((?<=[{ss_letters_crick}]))"
        f"([{ds_letters}]{{1,{length}}})"
        f"((?=[^{ss_letters_watson}{ds_letters}])))|"
        f"(?P<crick>((?<=[{ss_letters_watson}]))"
        f"([{ds_letters}]{{1,{length}}})"
        f"((?=[^{ss_letters_crick}{ds_letters}])))"
    )

    return re.compile(regex.encode("ascii"))


def regex_ds_melt_factory(length: int) -> re.Pattern:
    """
    A regular expression for finding double-stranded regions flanked by single-stranded DNA
    that can be melted to shed multiple double stranded fragments.

    This function returns a regular expression that finds double-stranded regions
    (of length <= length) that are flanked by single-stranded regions on opposite
    sides in dscode format. These regions are useful to identify as potential melt
    sites, since melting them leads to separation into multiple double stranded fragments.

    The regular expression finds double stranded patches flanked by empty
    positions on opposite sides(see figure below). Melting of this kind of
    sites leads to separation into multiple double stranded fragments.

    ::
        aaaGFTTAIAttt   <-- dscode

        aaaG TTACAttt   <-- "TTA" is found by the regex for length <= 3
        tttCTAAT Taaa

    Examples
    --------

    >>> from pydna.dseq import Dseq
    >>> regex = regex_ds_melt_factory(3)
    >>> s = Dseq("aaaGFTTAIAttt")
    >>> s
    Dseq(-13)
    aaaG TTACAttt
    tttCTAAT Taaa
    >>> mobj = regex.search(s._data)
    >>> mobj.groupdict()
    {'watson': None, 'crick': b'TTA'}

    Parameters
    ----------
    length : int
        Max length of double stranded region flanked by single stranded
        regions.

    Returns
    -------
    TYPE
        regular expression object.

    """

    regex = (
        f"(?P<watson>((?<=[{ss_letters_watson}])|^)"
        f"([{ds_letters}]{{1,{length}}})"
        f"((?=[^{ss_letters_watson}{ds_letters}])|$))|"
        f"(?P<crick>((?<=[{ss_letters_crick}])|^)"
        f"([{ds_letters}]{{1,{length}}})"
        f"((?=[^{ss_letters_crick}{ds_letters}])|$))"
    )

    return re.compile(regex.encode("ascii"))


def anneal_strands(strand_a: str, strand_b: str) -> bool:
    """
    Test if two DNA strands containing dscode anneal or not.

    Both strands are assumed to be given in 5' -> 3' direction.

    Examples
    --------

    >>> from pydna.alphabet import anneal_strands
    >>> a = "TTA"
    >>> b = "AAT"[::-1]
    >>> anneal_strands(a, b)
    True
    >>> anneal_strands(b, a)
    True
    >>> c = "UUA"
    >>> anneal_strands(c, b)
    True
    >>> anneal_strands(a.lower(), b)
    True
    >>> anneal_strands("TG", "AA")
    False

    Parameters
    ----------
    watson : str
        A single DNA strand.
    crick : str
        A single DNA strand.

    Returns
    -------
    bool
        True if annealing is perfect.

    """
    w = strand_a.translate(dscode_to_watson_table)
    c = strand_b.translate(complement_table_for_dscode).translate(
        dscode_to_crick_table
    )[::-1]
    for x, y in zip(w, c):
        try:
            basepair_dict[(x, y)]
        except KeyError:
            return False
    return True
