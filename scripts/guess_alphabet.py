#!/usr/bin/env python
# -*- coding: utf-8 -*-
def guess_alphabet(sequence):
    """
    guess_alphabet(sequence)

    returns an apropriate Biopython alphabet based on the content of the sequence.

    >>> from pydna.dsdna import guess_alphabet
    >>> guess_alphabet("aaa")
    IUPACUnambiguousDNA()
    >>> guess_alphabet("bbb")
    IUPACProtein()
    >>> guess_alphabet("aaaac")
    IUPACUnambiguousDNA()
    >>> guess_alphabet("agtc")
    IUPACUnambiguousDNA()
    >>> guess_alphabet("aguc")
    IUPACUnambiguousRNA()
    >>> guess_alphabet("ALAASN")
    'three letter code'

    """
    # SeqRecord or Dseqrecord?
    if hasattr(sequence, "features"):
        sequence = str(sequence.seq)

    # Seq or Dseq ?
    elif hasattr(sequence, "reverse_complement"):
        sequence = str(sequence)

    else:
        if not isinstance(sequence, str):
            warnings.warn(
                ("Input is {}, expected string, " "unicode string, Seq or SeqRecord " "objects!").format(
                    type(sequence)
                )
            )
        sequence = str(sequence)

    if len(sequence) < 1:
        return SingleLetterAlphabet()

    for c in sequence:
        if c not in string.printable:
            return SingleLetterAlphabet()

    xp = set(extended_protein.letters)
    pr = set(protein.letters)

    ad = set(ambiguous_dna.letters)
    ud = set(unambiguous_dna.letters)
    ed = set(extended_dna.letters)

    ar = set(ambiguous_rna.letters)
    ur = set(unambiguous_rna.letters)

    all = xp | pr | ad | ud | ed | ar | ur

    sequence_chars = set(sequence.upper())

    if sequence_chars - all - set(string.punctuation + string.whitespace):
        return SingleLetterAlphabet()

    nucleic_count = 0

    for letter in "GATCUNgatcun":
        nucleic_count += sequence.count(letter)

    if float(nucleic_count) / float(len(sequence)) >= 0.9:  # DNA or RNA
        if "T" in sequence_chars and "U" in sequence_chars:
            alphabet = NucleotideAlphabet()
        elif not sequence_chars - ud:
            alphabet = unambiguous_dna
        elif not sequence_chars - ad:
            alphabet = ambiguous_dna
        elif not sequence_chars - ed:
            alphabet = extended_dna
        elif not sequence_chars - ur:
            alphabet = unambiguous_rna
        elif not sequence_chars - ar:
            alphabet = extended_rna
        else:
            alphabet = NucleotideAlphabet()
    else:
        threecode = [
            "ALA",
            "ASX",
            "CYS",
            "ASP",
            "GLU",
            "PHE",
            "GLY",
            "HIS",
            "ILE",
            "LYS",
            "LEU",
            "MET",
            "ASN",
            "PRO",
            "GLN",
            "ARG",
            "SER",
            "THR",
            "VAL",
            "TRP",
            "TYR",
            "GLX",
            "XAA",
            "TER",
            "SEL",
            "PYL",
            "XLE",
        ]
        tc = set(threecode)
        three_letter_alphabet = set([sequence[i : i + 3] for i in range(0, len(sequence), 3)])
        if not three_letter_alphabet - tc:
            alphabet = "three letter code"
        elif sequence_chars - pr:
            alphabet = protein
        elif sequence_chars - xp:
            alphabet = extended_protein
        else:
            alphabet = ProteinAlphabet()
    return alphabet
