# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import cre_lox_integration
from pydna.cre_lox import LOXP_SEQUENCE

a = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
b = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)

res = cre_lox_integration(a, b)


print(res[0].seq)
