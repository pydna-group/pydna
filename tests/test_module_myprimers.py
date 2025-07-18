#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

monkeypatch = pytest.MonkeyPatch()


def test_PrimerList_init(monkeypatch, capsys):
    monkeypatch.setenv("pydna_primers", "primers_linux_line_endings.txt")

    from pydna.parsers import parse_primers

    primer_source = parse_primers("primers_linux_line_endings.txt")[::-1]

    from pydna import myprimers

    from importlib import reload

    reload(myprimers)

    pl1 = myprimers.PrimerList()

    assert pl1 == primer_source

    assert len(pl1) == 4

    from pydna.parsers import parse_primers
    from pydna import myprimers

    primer_source = parse_primers("primers_linux_line_endings.txt")[::-1]

    pl2 = myprimers.PrimerList(primer_source)

    pl3 = myprimers.PrimerList(path="primers_linux_line_endings.txt")

    assert len(pl2) == len(pl3) == 4

    assert pl1 == pl2 == pl3

    newlist = parse_primers(
        """
                            >abc
                            aaa
                            >efg
                            ttt
                            """
    )

    np = pl2.assign_numbers(newlist)

    assert [s.name for s in parse_primers(np)] == ["5_abc", "4_efg"]

    newlist = parse_primers(
        """
                            >abc
                            aaa
                            >efg
                            tttttttt
                            """
    )

    np = pl2.assign_numbers(newlist)

    assert [s.name for s in parse_primers(np)] == ["4_abc", "0_primer"]

    with pytest.raises(ValueError):
        # ValueError: efg 8-mer:5'-tttttttt-3' is not in list
        pl2.pydna_code_from_list(newlist)

    captured = capsys.readouterr()

    assert captured.out == ">abc 3-mer\naaa\n"

    import textwrap

    code = textwrap.dedent(
        """\
    from pydna.parsers import parse_primers

    p = {}

    p[0], p[1], p[2], p[3] = parse_primers('''

    >0_primer 8-mer
    tttttttt
    >1_primer 8-mer
    gggggggg
    >2_primer 8-mer
    cccccccc
    >3_primer 8-mer
    aaaaaaaa
    ''')"""
    )

    assert pl1.pydna_code_from_list(pl1) == code

    pl4 = myprimers.PrimerList(path="primers_linux_line_endings.txt")

    assert pl4.accessed_indices == []
    pl4[1] = primer_source[1]
    assert pl4.accessed == primer_source[1:2]

    assert pl4.accessed_indices == [1]

    pl4.accessed_indices == []
    assert pl4[2:4] == primer_source[2:4]
    assert pl4.accessed_indices == [1, 2, 3]

    with pytest.raises(ValueError):
        myprimers.PrimerList(identifier="/")

    with pytest.raises(ValueError):
        pl3[2] = primer_source[1]

    with pytest.raises(IndexError):
        pl3[999] = primer_source[1]

    with pytest.raises(ValueError):
        pl2.open_folder()

    from unittest import mock

    subp = mock.MagicMock()

    monkeypatch.setattr("sys.platform", "linux")
    monkeypatch.setattr("subprocess.run", subp)

    pl4.open_folder()

    subp.assert_called_with(["xdg-open", pl3.path.parent])


def test_check_primer_numbers(monkeypatch):
    from pydna import myprimers

    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.check_primer_numbers(pl) == [pl[5]]


def test_undefined_sequence(monkeypatch):
    from pydna import myprimers

    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.undefined_sequence(pl) == [pl[2]]


def test_find_duplicate_primers(monkeypatch):
    from pydna import myprimers

    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.find_duplicate_primers(pl) == [[pl[1], pl[3]]]


if __name__ == "__main__":
    pytest.main([__file__, "-vvv", "-s"])
