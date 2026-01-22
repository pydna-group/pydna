# -*- coding: utf-8 -*-

from pydna.readers import read


def test_longername_gives_deseq():
    myseq = read(open("broken_genbank_files/pGreenLantern1.gb", "r").read())
    assert myseq.circular
    myseq = read(open("broken_genbank_files/fakeGenBankFile.gb", "r").read())
    assert myseq.circular
