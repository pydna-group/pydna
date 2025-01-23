# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord

dsr = Dseqrecord("ATGCAAACAGTAATGATATAAT")
dsr.reverse_complement()
