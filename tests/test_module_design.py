#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Restriction import BamHI
import pytest

from pydna.amplify import pcr
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.design import circular_assembly_fragments
from pydna.parsers import parse
from pydna.primer import Primer

frags = parse(
    """
>49
ccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatcc
>50
ctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatc
>51
gtagtgaaacatacacgttgctcgggttcaccccggtccgttctgagtcga
"""
)


bam = Dseqrecord(BamHI.site)


def test_primer_design_all_Dseqrecords():
    with pytest.raises(ValueError):
        assembly_fragments(frags, 20)


def test_primer_design_all_pcr_products():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments(x, 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2]).seq


def test_primer_design_first_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([frags[0], x[1], x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2]).seq


def test_primer_design_second_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], frags[1], x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2]).seq


def test_primer_design_third_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], x[1], frags[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2]).seq


def test_primer_design_first_and_third_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([frags[0], x[1], frags[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2]).seq


def test_primer_design_same_first_and_third_Dseqrecord():
    from pydna.dseqrecord import Dseqrecord

    x = [primer_design(f) for f in frags]
    y = assembly_fragments([frags[0], x[1], frags[0]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_circular()[0]
    assert result.seguid() == (frags[0] + frags[1]).looped().seguid()

    a = Dseqrecord("ccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatcc", name="a")
    b = Dseqrecord(
        "ctgtcgagaaacttgtatccctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatcccaaggacacaatcgagctc",
        name="b",
    )

    z = Assembly((a, b, a), limit=20)

    result = z.assemble_linear()[0]
    assert (
        str(result.seq)
        == "ccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatccctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatcccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatcc"
    )
    result = z.assemble_circular()[0]
    assert result.seguid() == (frags[0] + frags[1]).looped().seguid()


def test_primer_design_linker_first():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([bam, x[0], x[1], x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (bam + frags[0] + frags[1] + frags[2]).seq


def test_primer_design_linker_second():  #
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], bam, x[1], x[2]], 20, 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq


def test_primer_design_linker_third():  #
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], x[1], bam, x[2]], maxlink=6, overlap=20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + bam + frags[2]).seq


def test_primer_design_linker_last():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], x[1], x[2], bam], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2] + bam).seq


def test_primer_design_linker_second_after_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([frags[0], bam, x[1], x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq


def test_primer_design_linker_second_before_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], bam, frags[1], x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq


def test_primer_design_linker_third_after_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], frags[1], bam, x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + bam + frags[2]).seq


def test_primer_design_two_fragments():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], x[1]], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1]).seq


def test_primer_design_four_fragments():
    x = [primer_design(f) for f in frags]
    fourth = Dseqrecord("TAAAAATAAAATTGTTGACAGCAGAAGTGATATAGAAATTTGTTAATTATTA")
    y = assembly_fragments(x + [fourth], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + frags[1] + frags[2] + fourth).seq


def test_primer_design_two_fragments_linker_in_between():  #
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0], bam, x[1]], 20, 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (frags[0] + bam + frags[1]).seq


def test_primer_design_two_fragments_flanking_linkers():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([bam, x[0], x[1], bam], 20)
    z = Assembly(y, limit=20)
    result = z.assemble_linear()[0]
    assert result.seq == (bam + frags[0] + frags[1] + bam).seq


def test_primer_design_one_fragment_flanking_linkers():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([bam, x[0], bam], 20)
    assert y[0].seq == (bam + frags[0] + bam).seq


def test_primer_Design():
    """test_primer_design"""

    a = Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    b = Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    c = Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")

    frags = assembly_fragments([primer_design(r) for r in (a, b, c)])

    asm = Assembly(frags)

    assert asm.assemble_linear()[0].seguid() == "ldseguid=WEE7xDdKIPkrEBPcSDsuAJIJVnY"

    frags = assembly_fragments([primer_design(r) for r in (a, b, c, a)])

    a2 = pcr(frags[-1].forward_primer, frags[0].reverse_primer, a)

    asm = Assembly((a2, frags[1], frags[2]))

    assert asm.assemble_circular()[0].seguid() == "cdseguid=85t6tfcvWav0wnXEIb-lkUtrl4s"


def test_primer_Design_with_linker():
    """test_primer_design"""

    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    linker = Dseqrecord("AAATTTCCCGGG")
    c = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    frags = assembly_fragments((primer_design(b), linker, primer_design(c)))

    asm1 = Assembly(frags)

    assert asm1.assemble_linear()[0].seguid(), (b + linker + c).seguid() == "l95igKB8iKAKrvvqE9CYksyNx40"

    frags = assembly_fragments((primer_design(b), linker, primer_design(c), primer_design(b)))

    b2 = pcr(frags[-1].forward_primer, frags[0].reverse_primer, b)

    asm2 = Assembly((b2, frags[1], frags[2]))

    assert (b + linker + c).looped().seguid() == asm2.assemble_circular()[0].seguid()

    assert (b + linker + c).looped().seguid() == "cdseguid=LqQ1_uMp2AmEZ_L2I1_njIMkVDc"


def test_primer_Design_given_fw_primer():
    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    a = primer_design(b, fp=Primer("agctactgactattag"))
    assert str(a.reverse_primer.seq) == "tagatcagtacagtc"


def test_primer_Design_given_rv_primer():
    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    a = primer_design(b, rp=Primer("tagatcagtacagtca"))
    assert str(a.forward_primer.seq) == "agctactgactattag"


def test_primer_Design_given_wrong_fw_primer():
    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    with pytest.raises(ValueError):
        primer_design(b, fp=Primer("agctactgactattagC"))


def test_primer_Design_given_wrong_rv_primer():
    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    with pytest.raises(ValueError):
        primer_design(b, rp=Primer("tagatcagtacagtcaC"))


def test_primer_Design_given_both_primers():
    b = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    with pytest.raises(ValueError):
        primer_design(b, fp=Primer("agctactgactattag"), rp=Primer("tagatcagtacagtca"))


def test_primer_Design_multiple_products():
    b = Dseqrecord("agctactgactattaggggttaagctactgactattaggggtttctgatcatctgatctactatctgactgtactgatcta")
    from pydna import _PydnaWarning

    with pytest.warns(_PydnaWarning):
        primer_design(b)


def test_circular_assembly_fragments():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments(x, 20, circular=True)
    z = Assembly(y, limit=20)
    result = z.assemble_circular()[0]
    assert (
        str(result.seq)
        == "tctgagtcgaccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatccctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatcgtagtgaaacatacacgttgctcgggttcaccccggtccgt"
    )


def test_circular_assembly_fragments2():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments((frags[0], x[1], x[2]), 20, circular=True)
    z = Assembly(y, limit=20)
    result = z.assemble_circular()[0]
    assert (
        str(result.seq)
        == "ccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatccctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatcgtagtgaaacatacacgttgctcgggttcaccccggtccgttctgagtcga"
    )


def test_too_short_template():
    from pydna.amplicon import Amplicon
    from pydna.design import primer_design
    from pydna.dseqrecord import Dseqrecord

    fragment = Dseqrecord("GCCACCATGG")
    assert primer_design(fragment) == Amplicon("")


def test_get_tm_and_primer_from_above_and_below():
    """Starting from higher and lower lengths gives the same value"""
    from pydna.design import _design_primer
    from pydna.tm import tm_default

    for f in frags:
        for temp in [40, 45, 50, 55]:
            assert _design_primer(temp, f, 13, tm_default) == _design_primer(temp, f, 13, tm_default, 40)


def test_use_estimate_function():
    from pydna.design import primer_design
    from pydna.tm import tm_default

    def tm_alt_lower(*args):
        return tm_default(*args) - 10

    def tm_alt_upper(*args):
        return tm_default(*args) + 10

    for f in frags:
        for temp in [40, 45, 50, 55]:
            amp_lower_with_estimate = primer_design(
                f, target_tm=temp, tm_func=tm_alt_lower, estimate_function=tm_default
            )
            amp_lower_no_estimate = primer_design(f, target_tm=temp, tm_func=tm_alt_lower)

            assert amp_lower_with_estimate == amp_lower_no_estimate

            amp_upper_with_estimate = primer_design(
                f, target_tm=temp, tm_func=tm_alt_upper, estimate_function=tm_default
            )
            amp_upper_no_estimate = primer_design(f, target_tm=temp, tm_func=tm_alt_upper)

            assert amp_upper_with_estimate == amp_upper_no_estimate


def test_primer_design_correct_value():
    from pydna.tm import tm_default

    for original_target_tm in range(60, 65):
        for frag in frags:
            amp = primer_design(frag, target_tm=original_target_tm, limit=15)
            fwd, rvs = amp.primers()
            possible_fwd = [Primer(frag[0:i].seq) for i in range(15, 40)]
            possible_rvs = [Primer(frag[-i:].reverse_complement().seq) for i in range(15, 40)]

            # Finds the closest forward primer
            fwd_diff = [abs(tm_default(f.seq) - original_target_tm) for f in possible_fwd]
            correct_fwd = possible_fwd[fwd_diff.index(min(fwd_diff))]
            assert str(fwd.seq) == str(correct_fwd.seq)

            # Uses that primers Tm as the target Tm for the reverse primer
            new_target_tm = tm_default(correct_fwd.seq)
            rvs_diff = [abs(tm_default(f.seq) - new_target_tm) for f in possible_rvs]
            correct_rvs = possible_rvs[rvs_diff.index(min(rvs_diff))]
            assert str(rvs.seq) == str(correct_rvs.seq)


if __name__ == "__main__":
    pytest.cmdline.main([__file__, "-v", "-s"])
