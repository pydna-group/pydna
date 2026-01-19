# -*- coding: utf-8 -*-
import pytest
from pydna import _PydnaWarning


def test_add_feature():
    from pydna.seq import Seq
    from Bio.Seq import Seq as BSeq
    from Bio.SeqRecord import SeqRecord as BSeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.seqrecord import SeqRecord

    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(3, 9)
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq="GGATCC")
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Seq("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Dseq("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=SeqRecord("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=BSeqRecord(BSeq("GGATCC")))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Dseqrecord("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    with pytest.raises(TypeError):
        s.add_feature(seq=Dseqrecord("GGGGGG"))
    s = SeqRecord("tttATGaaaTAAggg")
    s.add_feature(3, 12)
    assert s.features[0].qualifiers["label"] == ["orf9"]

    from Bio.Seq import Seq

    from pydna.seqrecord import SeqRecord

    a = SeqRecord(Seq("atgtaa"))

    a.add_feature(2, 4)

    assert (
        a.list_features().get_string()
        == "+-----+---------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n+-----+---------------+-----+-----+-----+-----+------+------+\n|   0 | L:ft2         | --> | 2   | 4   |   2 | misc |  no  |\n+-----+---------------+-----+-----+-----+-----+------+------+"
    )
    a.features[0].qualifiers
    del a.features[0].qualifiers["label"]
    assert (
        a.list_features().get_string()
        == "+-----+---------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n+-----+---------------+-----+-----+-----+-----+------+------+\n|   0 | nd            | --> | 2   | 4   |   2 | misc |  no  |\n+-----+---------------+-----+-----+-----+-----+------+------+"
    )
    a.features[0].qualifiers["note"] = ["AwesomeFeature"]
    assert (
        a.list_features().get_string()
        == "+-----+------------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note    | Dir | Sta | End | Len | type | orf? |\n+-----+------------------+-----+-----+-----+-----+------+------+\n|   0 | N:AwesomeFeature | --> | 2   | 4   |   2 | misc |  no  |\n+-----+------------------+-----+-----+-----+-----+------+------+"
    )


def test_stamp():
    from pydna.seqrecord import SeqRecord
    from Bio.Seq import Seq as bpSeq

    a = SeqRecord("attt")
    assert a.stamp() == "lsseguid=ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "pydna lsseguid=ot6JPLeAeMmfztW1736Kc6DAqlo" in a.annotations["comment"]

    a = SeqRecord("attt")
    a.annotations["comment"] = "something"
    assert a.stamp() == "lsseguid=ot6JPLeAeMmfztW1736Kc6DAqlo"

    assert "something" in a.annotations["comment"]

    a = SeqRecord("attt")
    assert a.stamp() == "lsseguid=ot6JPLeAeMmfztW1736Kc6DAqlo"
    a.seq = a.seq + "a"

    with pytest.warns(_PydnaWarning):
        a.stamp()

    assert "pydna lsseguid=ot6JPLeAeMmfztW1736Kc6DAqlo" in a.annotations["comment"]
    assert "pydna lsseguid=sFjwltMUzuIRRD_zch5hCn_8RV8" in a.annotations["comment"]

    a = SeqRecord(bpSeq("attt"))
    with pytest.raises(AttributeError):
        a.stamp()


def test___hash__():
    from Bio.Seq import Seq
    from pydna.seqrecord import SeqRecord

    s = SeqRecord(Seq("GGATCC"))
    t = SeqRecord(Seq("GGATCC"))
    u = SeqRecord(Seq("GGATCc"))
    assert hash(s) == hash(t) != hash(u)

    assert s == t
    assert s != u
    assert s != "1"

    assert s < u
    assert not s > u

    with pytest.raises(TypeError):
        assert s > "3"
    with pytest.raises(TypeError):
        assert s < "3"


def test_lcs():
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as BSeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.seqrecord import SeqRecord

    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation, ExactPosition

    s = SeqRecord(Seq("GGATCC"))

    expected = SeqFeature()
    expected.__dict__ = {
        "location": FeatureLocation(ExactPosition(0), ExactPosition(6), strand=1),
        "type": "read",
        "id": "<unknown id>",
        "qualifiers": {
            "label": ["sequence"],
            "ApEinfo_fwdcolor": ["#DAFFCF"],
            "ApEinfo_revcolor": ["#DFFDFF"],
        },
    }
    assert s.lcs("GGATCC").__dict__ == {
        "location": None,
        "type": "",
        "id": "<unknown id>",
        "qualifiers": OrderedDict(),
    }
    assert s.lcs("GGATCC", limit=4).__dict__ == expected.__dict__
    assert s.lcs(Seq("GGATCC"), limit=4).__dict__ == expected.__dict__
    assert (
        s.lcs(BSeqRecord(Seq("GGATCC"), name="sequence"), limit=4).__dict__
        == expected.__dict__
    )
    assert s.lcs(Dseq("GGATCC"), limit=4).__dict__ == expected.__dict__
    assert (
        s.lcs(Dseqrecord(Dseq("GGATCC"), name="sequence"), limit=4).__dict__
        == expected.__dict__
    )
    assert (
        s.lcs(Dseqrecord("GGATCC", name="sequence"), limit=4).__dict__
        == expected.__dict__
    )


def test_format():
    from Bio.Seq import Seq
    from pydna.seqrecord import SeqRecord

    s = SeqRecord(Seq("GGATCC"))
    s.format("gb")
    s.format("genbank")
    s.format("fasta")


def test_seqrecord():
    from pydna import seqrecord
    from pydna.seq import Seq

    s = seqrecord.SeqRecord("gatt")
    assert s.name == s.locus == "name"
    assert s.id == s.accession == "id"
    assert s.definition == s.description == "description"

    t = s[2:]

    assert t.name == t.locus == "part_name"
    assert t.id == t.accession == "part_id"
    assert t.definition == t.description == "description"

    s = seqrecord.SeqRecord("gatt", name="1234567890123456")
    assert s.name == s.locus == "1234567890123456"
    assert s.id == s.accession == "id"
    assert s.definition == s.description == "description"

    t = s[2:]

    assert t.name == t.locus == "part_12345678901"
    assert t.id == t.accession == "part_id"
    assert t.definition == t.description == "description"

    s = Seq("ATGAAATAA")

    obj = seqrecord.SeqRecord(s, name="1234567890123456")
    assert obj.name == "1234567890123456"

    # with pytest.warns(_PydnaWarning, match="truncated") as pdw:
    #     obj = seqrecord.SeqRecord(s, name="12345678901234567")

    obj = seqrecord.SeqRecord(s, annotations={"date": "24-DEC-1970"})
    assert obj.annotations["date"] == "24-DEC-1970"

    obj = seqrecord.SeqRecord(s)

    assert obj.map_target is None

    assert obj.isorf()

    not_orf = Seq("TGAAATAAA")

    obj = seqrecord.SeqRecord(not_orf)

    assert not obj.isorf()

    s3 = Seq("aaaATGAAATAAttt")

    from Bio.SeqFeature import SeqFeature, FeatureLocation

    loc = FeatureLocation(3, 12)

    sf = SeqFeature(loc, type="misc", qualifiers={"label": "lbl"})

    obj = seqrecord.SeqRecord(s3, features=[sf])

    obj.add_colors_to_features_for_ape()

    assert obj.features[0].qualifiers["ApEinfo_fwdcolor"] == ["#66ffa3"]
    assert obj.features[0].qualifiers["ApEinfo_revcolor"] == ["#66ffff"]

    assert obj.gc() == 0.067

    assert repr(obj) == (
        "SeqRecord(seq=Seq('aaaATGAAATAAttt'), id='id', "
        "name='name', description='description', dbxrefs=[])"
    )

    obj.annotations = {"date": "24-DEC-1970"}

    st = str(
        "ID: id\n"
        "Name: name\n"
        "Description: description\n"
        "Number of features: 1\n"
        "/date=24-DEC-1970\n"
        "Seq('aaaATGAAATAAttt')"
    )

    assert str(obj) == st

    obj.locus = "new123"
    assert obj.locus == obj.name == "new123"
    obj.locus = "new456"
    assert obj.locus == obj.name == "new456"

    obj.accession = "new123"
    assert obj.accession == obj.id == "new123"
    obj.id = "new456"
    assert obj.accession == obj.id == "new456"

    obj.definition = "new123"
    assert obj.definition == obj.description == "new123"
    obj.description = "new456"
    assert obj.definition == obj.description == "new456"

    with pytest.warns(_PydnaWarning, match="truncated"):
        obj.locus = "12345678901234567"

    lf = str(
        "+-----+---------------+-----+-----+-----+-----+------+------+\n"
        "| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n"
        "+-----+---------------+-----+-----+-----+-----+------+------+\n"
        "|   0 | L:l b l       | --- | 3   | 12  |   9 | misc | yes  |\n"
        "+-----+---------------+-----+-----+-----+-----+------+------+"
    )

    assert obj.list_features().get_string() == lf

    exft = obj.extract_feature(0)

    assert str(exft.seq) == "ATGAAATAA"
    import textwrap

    textwrap.dedent(
        """LOCUS       1234567890123456          15 bp    DNA              UNK 24-DEC-1970
    DEFINITION  new456.
    ACCESSION   new456
    VERSION     new456
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
         misc            4..12
                         /label="lbl"
                         /ApEinfo_fwdcolor="#7f3f3f"
                         /ApEinfo_revcolor="#7f3f3f"
    ORIGIN
            1 aaaatgaaat aattt
    //""".strip()
    )

    # print()
    # print(textwrap.dedent(...))

    # print(str(obj.format("gb")))

    # assert gbf+"\n" == str(obj.format("gb"))

    # print(obj.__hash__())

    s = Seq("atgCGACGGCGCCCGCTCGCGtaa")

    obj = seqrecord.SeqRecord(s, name="1234567890123456")

    assert obj.name == "1234567890123456"

    assert seqrecord.SeqRecord("atgtaa").gc() == 0.167

    assert seqrecord.SeqRecord("atgtaa").startcodon() == 1.0

    assert seqrecord.SeqRecord("atgtaa").stopcodon() == 0.47

    assert seqrecord.SeqRecord("atgtaa").gc() == 0.167


def test_cai():
    pytest.importorskip("cai2")
    from pydna import seqrecord
    from pydna.seq import Seq
    from pydna.codon import rare_codons

    assert seqrecord.SeqRecord("atgtaa").cai() == 1.0

    organism = "sce"

    s = seqrecord.SeqRecord(Seq("atg" + "".join(rare_codons[organism]) + "taa"))

    features = s.rarecodons()

    for feat in features:
        assert feat.extract(s).seq in s

    lol = [
        [
            "cds",
            "len",
            "cai",
            "gc",
            "sta",
            "stp",
            "n-end",
            "CGA",
            "CGG",
            "CGC",
            "CCG",
            "CTC",
            "GCG",
            "rare",
        ],
        ["ATG...TAA", 8.0, 0.219, 0.708, 1.0, 0.47, "2 min", 1, 1, 1, 1, 1, 1, 0.75],
    ]

    assert s.express().lol() == lol


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])


#    >>> import warnings
#    >>> from pydna import _PydnaWarning
#    >>> warnings.simplefilter('ignore', _PydnaWarning)

# obj.name = "1111111111111111"


#    def __init__(self, seq, id="<unknown id>", name="<unknown name>",
#                 description="<unknown description>", dbxrefs=None,
#                 features=None, annotations=None,
#                 letter_annotations=None):
#        """Create a SeqRecord.
