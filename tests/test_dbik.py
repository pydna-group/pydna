#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature,
    CompoundLocation,
    SimpleLocation,
    ExactPosition,
)
from Bio.SeqRecord import SeqRecord
from pydna.assembly import Assembly
from pydna.assembly2 import Assembly as Assembly2
from pydna.dseqrecord import Dseqrecord


def test_cut_dbik():

    cds = (
        "ATGGCATGCAAATTACCGGTAAGCATAGCTAATTGCAGTAACCGTTGAATAGCATTACGT"
        "AGCTAGTACGATCAATCGGATGCATAGCATGCATAGCAGTACGTACGGCATAGCTAGCAT"
        "AGCAGCATCGTACGTACGTACGATGCATGCATGCATGCATTGGCAATGCATGCATGCATA"
    )
    body = (
        "TTAACGGTACATCGCATGCATGCATAGCATAGCATCGGTAGCATGCATAGCAGCATAGCA"
        "TGCATAGCATGCATGCATGCATAGCATAGCATAGCATCGGGCATGCATGCATGCATGCAT"
        "AGCATAGCATGCATGCATGCATGCAACTAGCATCGTACGTACGCATGCAGCATCGGGCAT"
    )
    ovl = 30

    # Backbone whose 5' end carries the CDS tail and 3' end carries the CDS head,
    # forcing the CDS to span the assembled origin once the insert is folded in.
    backbone = cds[-ovl:] + body + cds[:ovl]

    # Insert is the CDS in RC form. Assembly must reverse-complement it to match
    # overlaps.
    insert_rc = str(Seq(cds).reverse_complement())

    bb = SeqRecord(
        Seq(backbone),
        id="bb",
        name="bb",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    ins = SeqRecord(
        Seq(insert_rc),
        id="ins",
        name="ins",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    ins.features.append(
        SeqFeature(
            FeatureLocation(0, len(cds), strand=-1),
            type="CDS",
            qualifiers={"label": ["myCDS"]},
        )
    )

    asm = Assembly([Dseqrecord(bb), Dseqrecord(ins)], limit=ovl)
    asm2 = Assembly2([Dseqrecord(bb), Dseqrecord(ins)], limit=ovl)

    prod1, prod2 = asm.assemble_circular()
    prod3, *_ = asm2.assemble_circular()

    prod1.name = "bb_ins1"
    prod2.name = "bb_ins2"
    prod3.name = "bb_ins3"

    cds_feat1 = next(
        f for f in prod1.features if f.qualifiers.get("label", [""])[0] == "myCDS"
    )

    assert cds_feat1.location == CompoundLocation(
        [
            SimpleLocation(ExactPosition(210), ExactPosition(360), strand=1),
            SimpleLocation(ExactPosition(0), ExactPosition(30), strand=1),
        ],
        "join",
    )
    assert cds_feat1.location.strand == 1
    assert [p.strand for p in cds_feat1.location.parts] == [1, 1]

    cds_feat2 = next(
        f for f in prod2.features if f.qualifiers.get("label", [""])[0] == "myCDS"
    )

    assert cds_feat2.location == SimpleLocation(
        ExactPosition(0), ExactPosition(180), strand=-1
    )
    assert cds_feat2.location.strand == -1
    assert [p.strand for p in cds_feat2.location.parts] == [-1]

    cds_feat3 = next(
        f for f in prod3.features if f.qualifiers.get("label", [""])[0] == "myCDS"
    )

    assert cds_feat3.location == CompoundLocation(
        [
            SimpleLocation(ExactPosition(210), ExactPosition(360), strand=1),
            SimpleLocation(ExactPosition(0), ExactPosition(30), strand=1),
        ],
        "join",
    )
    assert cds_feat3.location.strand == 1
    assert [p.strand for p in cds_feat3.location.parts] == [1, 1]

    assert cds_feat3 == cds_feat1

    assert (
        prod1.seguid()
        == prod2.seguid()
        == prod3.seguid()
        == "cdseguid=EN8cEkWMLtmYLWqha4Kf2dQHz8A"
    )
