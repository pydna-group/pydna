# -*- coding: utf-8 -*-
#                  _
#  _ __  _   _  __| |_ __   __ _
# | '_ \| | | |/ _` | '_ \ / _` |
# | |_) | |_| | (_| | | | | (_| |
# | .__/ \__, |\__,_|_| |_|\__,_|
# |_|    |___/
#
#
#                                .:::.
#                     .^!???!::?55YYY55?^
#                   .7PY7!!!?#B?:    .:?GY7!~^:.
#                  .5G^    .?#~ :^.     ^#P?JY55Y7^
#                 :J&7   7GG&B J#BB!    .GY^^^~~7J557.
#               ~5PYG5.  !PG5BJ!Y5J:   .JG7^~~~^^^^!YG7
#             .5GJ!~!5P?^.:^^Y#5?!^^^!?P5!~~~~~~~~~^^7G?
#            :GP!!!!!!7YY55YYJ!!?YYYYY?!~~~~~~~~~~~~~^?B~
#            5G!!!!!!!!!!!!!~~!!~~~~~~~~~~~!?7~~~~~~~^!B7
#           ^#J!!!!!!!!!!!!!!!!!!~~~~~!!7?JJ7!~~~~~~~~YP:
#           :B577777777777777777???JJJJJ?7!~~~~~~~~~~5P:
#          ^JGPB#GP5YYYJJJJYYYJJJJJ?7!!!~~!!!~~~~~~!PP:
#        :5P7: :~7?JYYYY5555YYYYYYJ5B7!!!!!!!!!!!~7BY.
#       ~BY.           .......     !#?!!!!!!!!!!~JB7
#      !B?                         YB7!!!!!!!!!!5B~                        .:.
#  .^75&J                        .JB?!!!!!!!!!7G5.                     .^75#5:
# .7Y?~GP.                      :5G7!777!!!!!?BJ             .:^~~!7??Y55G&Y.
#      ~?.                     !G5777777777!JB7           :!Y55YYYYYJJ??PG!
#                            .JBJ777777777!J#G7!^:.    .!YP5?7!!!!!7?JPP?.
#                           ~GG?77777777775GJ?Y5555Y?!7P5?!!!!7?JYPP5J!.
#                      .^!?5GY7777777777JGP777!!!77?YP#P7!77?YPPY?~:
#                   :?Y55#&5?7???77777?5GY7777777777!!7YG5J5PJ~.
#                 ~5PY7?GP??????????7?55?77777777777777!75&&?.
#               :5BJ7!?B5?????????????777777777777777777!?BYYP?.
#              ^BP77775BJ????????????????????77777777777JG5~^~PP.
#             .GG????7?GGJ??????????????????????77777?YPP?^^~^^BY
#             ?#????????5GPYJ????????????????????JY5P5Y!^^^^^^^PG
#             JB??????????YPPPPP555YYYYY5555555555J?!~^^^^^^^^!BB?
#             7&P??????????????JJYYY55YYYJJ??7!!~~~~~~~~~~^^^?GJ~G5.
#            ^#YYG5J????????????77777777!!!!!!!!!!!~~~~~~~!J55!^^!#?
#            5B~~!YPP5YJJ??????????77777777!!!!!!!!!77?JY55J!^^^^^P5
#           .B5!!!!!7JYPPPPP5555555555555YYYY55Y55555YJ?!~^^^^^~^~#J
#            5B7!777!!!!!77???JJJJJJJJJJJJ???777!!!~~~~~~~~~~~~^!GP.
#            .5GJ7777777777777777!!!!!!!!!!!!!!!!!!!!!!!!!!~~!7YPJ.
#              !5P5Y??????????????77777777777777777!!!!!77?J5P5?:
#                ^7YPPP55YYJJJ???????????????JJJJJYY555555Y?~:
#                   .:^!?JY5555PPPPPPPPPP555555YYJ??7~^:.
#                            ...............


# Cloning with Python

import pydna
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import BamHI
from Bio.Restriction import EcoRI
from pydna.genbank import genbank
from pydna.parsers import parse_primers
from pydna.readers import read
from pydna.amplify import pcr
from pydna.assembly import Assembly

myinsertdna = "gtctgtgttgtt"

myplasmid = "ccatttgtatgttcagctaaCCCGGGcttctacccatcccccgaag"

SmaI = "CCCGGG"

myrecplasmid = myplasmid[:23] + myinsertdna + myplasmid[23:]

myrecplasmid = (
    myplasmid[: myplasmid.find(SmaI) + 3]
    + myinsertdna
    + myplasmid[myplasmid.find(SmaI) + 3 :]
)


pydna.logo()

Dseq("aaa", "ttt")

Dseq("aaa", "ttt", ovhg=0)

Dseq("aaa", "ttt", ovhg=1)

Dseq("aaa", "ttt", ovhg=-1)

Dseq("aaa", "ttt", linear=False, ovhg=0)

Dseq("aaa", "ttt", circular=True, ovhg=0)

Dseq("ggaaa", "ccttt")


x = Dseqrecord("GGATCC")

x.seq

a, b = Dseqrecord("GGATCC").cut(BamHI)

a, b

a.seq

b.seq

(a + b).seq

(b + a).seq

c, d = Dseqrecord("GAATTC").cut(EcoRI)

c.seq

d.seq

gbfile = genbank("CS570233.1")

gbfile

a, b = parse_primers(
    """
>myprimer1
gtcatctacgtcgtacgt
>myprimer2
gtgtaggtctatttagtcgtag
"""
)

mytemplate = read(
    """
>mytemplate
gtcatctacgtcgtacgttgtgtgtacgtagtagtgtcactacgactaaatagacctacac
"""
)

ampl = pcr(a, b, mytemplate)

a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
x = Assembly((a, b, c), limit=14)

candidates = x.assemble_circular()

candidates

x, y = candidates

x
