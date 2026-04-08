# -*- coding: utf-8 -*-

from Bio import SeqIO
from io import StringIO
from pydna.parsers import parse_proteins
from pydna.align import align

fasta_sequence = """\
>sp|P12345|MY_PROT Some protein description
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQ
"""

genpept_sequence = """\
LOCUS       AF171097_1               143 aa            linear   BCT 21-AUG-2001
DEFINITION  transcriptional regulator RovA [Yersinia enterocolitica].
ACCESSION   AAD51968
VERSION     AAD51968.1
DBSOURCE    accession AF171097.1
KEYWORDS    .
SOURCE      Yersinia enterocolitica
  ORGANISM  Yersinia enterocolitica
            Bacteria; Pseudomonadati; Pseudomonadota; Gammaproteobacteria;
            Enterobacterales; Yersiniaceae; Yersinia.
REFERENCE   1  (residues 1 to 143)
  AUTHORS   Revell,P.A. and Miller,V.L.
  TITLE     A chromosomally encoded regulator is required for expression of the
            Yersinia enterocolitica inv gene and for virulence
  JOURNAL   Mol. Microbiol. 35 (3), 677-685 (2000)
   PUBMED   10672189
REFERENCE   2  (residues 1 to 143)
  AUTHORS   Revell,P.A. and Miller,V.L.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-JUL-1999) Molecular Microbiology, Washington
            University School of Medicine, Campus Box 8230, 660 South Euclid,
            St. Louis, MO 63110, USA
COMMENT     Method: conceptual translation.
FEATURES             Location/Qualifiers
     source          1..143
                     /organism="Yersinia enterocolitica"
                     /strain="JB580v"
                     /serotype="O:8"
                     /db_xref="taxon:630"
     Protein         1..143
                     /product="transcriptional regulator RovA"
                     /name="regulates inv expression"
     Region          1..143
                     /region_name="PRK03573"
                     /note="transcriptional regulator SlyA; Provisional"
                     /db_xref="CDD:179596"
     CDS             1..143
                     /gene="rovA"
                     /coded_by="AF171097.1:380..811"
                     /note="regulator of virulence"
                     /transl_table=11
ORIGIN
        1 mestlgsdla rlvrvwrali dhrlkplelt qthwvtlhni nrlppeqsqi qlakaigieq
       61 pslvrtldql eekglitrht candrrakri klteqsspii eqvdgvicst rkeilggisp
      121 deiellsgli dklerniiql qsk
//
"""


def read_w_seqio(sequence: str, sequenceformat: str):
    handle = StringIO(sequence)
    record = next(SeqIO.parse(handle, sequenceformat))
    return record


def test_fasta_protein():
    (pydna_protein,) = parse_proteins(fasta_sequence)
    biopython_protein = read_w_seqio(fasta_sequence, "fasta-blast")
    assert pydna_protein.__dict__.keys() == biopython_protein.__dict__.keys() | {
        "map_target"
    }


def test_genpept_protein():
    (pydna_protein,) = parse_proteins(genpept_sequence)
    biopython_protein = read_w_seqio(genpept_sequence, "gb")
    assert pydna_protein.__dict__.keys() == biopython_protein.__dict__.keys() | {
        "map_target"
    }


def test_same_protein():

    target = "MKTAYIAKKKKKISFVKSHFSR"
    query = "MKTAYIAKKKKKISFVKSHFSR"

    aln, editlist = align(target, query)

    assert str(aln) == (
        "target            0 MKTAYIAKKKKKISFVKSHFSR 22\n"
        "                  0 |||||||||||||||||||||| 22\n"
        "query             0 MKTAYIAKKKKKISFVKSHFSR 22\n"
    )

    assert editlist == []


def test_protein_inserts():

    target = "MKTAYIAKQRQISFVKSHFSRQ"
    query = "MKTAYIAKQISFVKSHFSR"

    aln, editlist = align(target, query)

    assert str(aln) == (
        "target            0 MKTAYIAKQRQISFVKSHFSRQ 22\n"
        "                  0 ||||||||--|||||||||||- 22\n"
        "query             0 MKTAYIAK--QISFVKSHFSR- 19\n"
    )

    assert editlist == ["Insert QR at position 9-10", "Insert Q at position 22"]


def test_protein3():

    target = "MKTAYIAKQISFVKSHFSR"
    query = "MKTAYIAKQRQISFVKSHFSRQ"

    aln, editlist = align(target, query)

    assert str(aln) == (
        "target            0 MKTAYIAK--QISFVKSHFSR- 19\n"
        "                  0 ||||||||--|||||||||||- 22\n"
        "query             0 MKTAYIAKQRQISFVKSHFSRQ 22\n"
    )

    assert editlist == ["Delete QR after position 8", "Delete Q after position 19"]


def test_protein_insert_substitute_delete():

    target = "MKTAYIAKKKKKISFVKSHFSR"
    query = "MKTAYIAKQRQISFVKSHFSRQ"

    aln, editlist = align(target, query)

    assert str(aln) == (
        "target            0 MKTAYIAKKKKKISFVKSHFSR- 22\n"
        "                  0 |||||||-|...||||||||||- 23\n"
        "query             0 MKTAYIA-KQRQISFVKSHFSRQ 22\n"
    )

    assert editlist == [
        "Insert K at position 8",
        "Substitute QRQ → KKK from position 10 to 12",
        "Delete Q after position 22",
    ]
