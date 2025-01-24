# -*- coding: utf-8 -*-
"""docstring.

https://en.wikipedia.org/wiki/SHA-1

    raw data      : The quick brown fox jumps over the lazy dog
    SHA1 (dec)    : 273069992013452546326057769888623105462687230738
    SHA1 (hex)    : 2fd4e1c67a2d28fced849ee1bb76e7391b93eb12
    SHA1 (base64) : L9ThxnotKPzthJ7hu3bnORuT6xI

"""

import hashlib
import base64

text = "The quick brown fox jumps over the lazy dog"

SHA1_hex = hashlib.sha1(text.encode("ASCII")).hexdigest()

assert SHA1_hex == "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12"

SHA1_dec = int(SHA1_hex, 16)

SHA1_bin = bin(SHA1_dec)[2:]

bytestring = hashlib.sha1(text.encode("ASCII")).digest()

b64 = base64.standard_b64encode(bytestring).decode("ASCII")

print("raw data      :", text)
print("SHA1 (bin)    :", SHA1_bin[:10] + " .. " + SHA1_bin[-10:] + f" ({len(SHA1_bin)} bits)")
print("SHA1 (dec)    :", SHA1_dec)
print("SHA1 (hex)    :", SHA1_hex)
print("SHA1 (base64) :", b64.rstrip("="))


text = "The quick brown foxes jump over the lazy dog."

bytestring = hashlib.sha1(text.encode()).digest()

SHA1_hex = hashlib.sha1(text.encode("ASCII")).hexdigest()

SHA1_dec = int(SHA1_hex, 16)

SHA1_bin = bin(SHA1_dec)[2:]

b64 = base64.standard_b64encode(bytestring).decode("ASCII")

b64us = base64.urlsafe_b64encode(bytestring).decode("ASCII")

print()
print("raw data             :", text)
print("SHA1 (bin)           :", SHA1_bin[:10] + " .. " + SHA1_bin[-10:] + " (160 bits)")
print("SHA1 (dec)           :", SHA1_dec)
print("SHA1 (hex)           :", SHA1_hex)
print("SHA1 (base64)        :", b64.rstrip("="))
print("SHA1 (base64urlsafe) :", b64us.rstrip("="))

"""

raw data             : The quick brown foxes jump over the lazy dog.
SHA1 (dec)           : 1036398307098929586513334149440434155590523968773
SHA1 (hex)           : b589b517fdceaf32032fe5f4caf0805d536e5d05
SHA1 (base64)        : tYm1F/3OrzIDL+X0yvCAXVNuXQU
SHA1 (base64urlsafe) : tYm1F_3OrzIDL-X0yvCAXVNuXQU
"""


text = "GATTACAGATTACAGATTACA"

bytestring = hashlib.sha1(text.encode()).digest()

SHA1_hex = bytestring.hex()

SHA1_dec = int(SHA1_hex, 16)

SHA1_bin = bin(SHA1_dec)[2:]

b64 = base64.standard_b64encode(bytestring).decode("ASCII")

b64us = base64.urlsafe_b64encode(bytestring).decode("ASCII")

print()
print("raw data             :", text)
print("SHA1 (bin)           :", SHA1_bin[:10] + " .. " + SHA1_bin[-10:] + " (160 bits)")
print("SHA1 (dec)           :", SHA1_dec)
print("SHA1 (hex)           :", SHA1_hex)
print("SHA1 (base64)        :", b64.rstrip("="))
print("SHA1 (base64urlsafe) :", b64us.rstrip("="))


text = "GATTACAGATTACAGATTACA"

bytestring = hashlib.sha256(text.encode()).digest()

SHA256_hex = bytestring.hex()

SHA256_dec = int(SHA256_hex, 16)

SHA256_bin = bin(SHA256_dec)[2:]

b64 = base64.standard_b64encode(bytestring).decode("ASCII")

b64us = base64.urlsafe_b64encode(bytestring).decode("ASCII")

print()
print("raw data               :", text)
print("SHA256 (bin)           :", SHA256_bin[:10] + " .. " + SHA256_bin[-10:] + f" ({len(SHA256_bin)} bits)")
print("SHA256 (dec)           :", SHA256_dec)
print("SHA256 (hex)           :", SHA256_hex)
print("SHA256 (base64)        :", b64.rstrip("="))
print("SHA256 (base64urlsafe) :", b64us.rstrip("="))


"""

dsDNA  ovhg    rc        ovhg    raw               SEGUID

  NNN   2      NNNNN     0        0NNNNN|NNN       Pi78amLZSO-FUb8yitZI6kaEXEs
NNNNN          NNN

 NNNN   1      NNNNN     0        0NNNNN|NNNN      Q17Z6aJWSFjK269FiQ17dodqRVY
NNNNN          NNNN

NNNNN  -1      NNNN      0       -1NNNNN|NNNN      ODxYaIo5ZF65_0e23SQgalXHrMA
 NNNN          NNNNN

NNNNN  -2      NNN       0       -2NNNNN|NNN       H6qxFiotG_kLdpZDwGJkqc2tExI
  NNN          NNNNN

"""


from pydna.dseq import Dseq

frags = (("NNN", "NNNNN", 2), ("NNNN", "NNNNN", 1), ("NNNNN", "NNNNN", 0), ("NNNNN", "NNNN", -1), ("NNNNN", "NNN", -2))

for frag in frags:
    x = Dseq(*frag)
    print(repr(x))
    print(x.lseguid())
