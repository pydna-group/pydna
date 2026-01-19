#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""This module provide functions for melting temperature calculations."""


import math
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

import textwrap
from pydna._pretty import pretty_str as ps

# See the documentation for Bio.SeqUtils.MeltingTemp for more details
# The 10X Taq Buffer with (NH4)2SO4 is commercialized by companies like
# ThermoFisher, although we make it ourselves
# 10X Buffer Composition
# 750 mM Tris-HCl (pH 8.8 at 25°C),
# 200 mM (NH4)2SO4,
# 0.1% (v/v) Tween 20.


def tm_default(
    seq,
    check=True,
    strict=True,
    c_seq=None,
    shift=0,
    nn_table=mt.DNA_NN4,  # DNA_NN4: values from SantaLucia & Hicks (2004)
    tmm_table=None,
    imm_table=None,
    de_table=None,
    dnac1=500 / 2,  # I assume 500 µM of each primer in the PCR mix
    dnac2=500 / 2,  # This is what MELTING and Primer3Plus do
    selfcomp=False,
    Na=40,
    K=0,
    Tris=75.0,  # We use the 10X Taq Buffer with (NH4)2SO4 (above)
    Mg=1.5,  # 1.5 mM Mg2+ is often seen in modern protocols
    dNTPs=0.8,  # I assume 200 µM of each dNTP
    saltcorr=7,  # Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log[Na+]
    func=mt.Tm_NN,  # Used by Primer3Plus to calculate the product Tm.
):
    return func(
        seq,
        check=check,
        strict=strict,
        c_seq=c_seq,
        shift=shift,
        nn_table=nn_table,
        tmm_table=tmm_table,
        imm_table=imm_table,
        de_table=de_table,
        dnac1=dnac1,
        dnac2=dnac2,
        selfcomp=selfcomp,
        Na=Na,
        K=K,
        Tris=Tris,
        Mg=Mg,
        dNTPs=dNTPs,
        saltcorr=saltcorr,
    )


def tm_dbd(
    seq,
    check=True,
    strict=True,
    c_seq=None,
    shift=0,
    nn_table=mt.DNA_NN3,
    tmm_table=None,
    imm_table=None,
    de_table=None,
    dnac1=250,
    dnac2=250,
    selfcomp=False,
    Na=50,
    K=0,
    Tris=0,
    Mg=1.5,
    dNTPs=0.8,
    saltcorr=1,
    func=mt.Tm_NN,
):
    return func(
        seq,
        check=check,
        strict=strict,
        c_seq=c_seq,
        shift=shift,
        nn_table=nn_table,
        tmm_table=tmm_table,
        imm_table=imm_table,
        de_table=de_table,
        dnac1=dnac1,
        dnac2=dnac2,
        selfcomp=selfcomp,
        Na=Na,
        K=K,
        Tris=Tris,
        Mg=Mg,
        dNTPs=dNTPs,
        saltcorr=saltcorr,
    )


def tm_product(seq: str, K=0.050):
    """Tm calculation for the amplicon.

    according to:

    Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
    ing temperature for DNA amplification in vitro
    http://www.ncbi.nlm.nih.gov/pubmed/2243783
    """
    tmp = 81.5 + 0.41 * gc_fraction(seq) * 100 + 16.6 * math.log10(K) - 675 / len(seq)
    return tmp


def ta_default(fp: str, rp: str, seq: str, tm=tm_default, tm_product=tm_product):
    """Ta calculation.

    according to:

    Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
    ing temperature for DNA amplification in vitro
    http://www.ncbi.nlm.nih.gov/pubmed/2243783

    The formula described uses the length and GC content of the product and
    salt concentration (monovalent cations)
    """
    return 0.3 * min((tm(fp), tm(rp))) + 0.7 * tm_product(seq) - 14.9


def ta_dbd(fp, rp, seq, tm=tm_dbd, tm_product=None):
    return min((tm(fp), tm(rp))) + 3


def program(amplicon, tm=tm_default, ta=ta_default):
    r"""Returns a string containing a text representation of a suggested
    PCR program using Taq or similar polymerase.

    ::

     |95°C|95°C               |    |tmf:59.5
     |____|_____          72°C|72°C|tmr:59.7
     |3min|30s  \ 59.1°C _____|____|60s/kb
     |    |      \______/ 0:32|5min|GC 51%
     |    |       30s         |    |1051bp

    """

    taq_extension_rate = 45  # seconds/kB PCR product length (1min/kb)
    extension_time_taq = max(30, int(taq_extension_rate * len(amplicon) / 1000))
    # seconds

    f = textwrap.dedent(
        r"""
        |95°C|95°C               |    |tmf:{tmf:.1f}
        |____|_____          72°C|72°C|tmr:{tmr:.1f}
        |3min|30s  \ {ta:.1f}°C _____|____|{rate}s/kb
        |    |      \______/{0:2}:{1:0>2}|5min|GC {GC}%
        |    |       30s         |    |{size}bp
        """.format(
            rate=taq_extension_rate,
            size=len(amplicon.seq),
            ta=round(
                ta(
                    amplicon.forward_primer.footprint,
                    amplicon.reverse_primer.footprint,
                    str(amplicon.seq),
                ),
                1,
            ),
            tmf=tm(amplicon.forward_primer.footprint),
            tmr=tm(amplicon.reverse_primer.footprint),
            GC=int(amplicon.gc() * 100),
            *map(int, divmod(extension_time_taq, 60)),
        )
    ).strip()

    return ps(f)


taq_program = program


def dbd_program(amplicon, tm=tm_dbd, ta=ta_dbd):
    r"""Text representation of a suggested PCR program.

    Using a polymerase with a DNA binding domain such as Pfu-Sso7d.

    ::

     |98°C|98°C               |    |tmf:53.8
     |____|_____          72°C|72°C|tmr:54.8
     |30s |10s  \ 57.0°C _____|____|15s/kb
     |    |      \______/ 0:15|5min|GC 51%
     |    |       10s         |    |1051bp


     |98°C|98°C      |    |tmf:82.5
     |____|____      |    |tmr:84.4
     |30s |10s \ 72°C|72°C|15s/kb
     |    |     \____|____|GC 52%
     |    |      3:45|5min|15058bp

    """
    PfuSso7d_extension_rate = 15  # seconds/kB PCR product
    extension_time_PfuSso7d = max(
        10, int(PfuSso7d_extension_rate * len(amplicon) / 1000)
    )  # seconds

    # The program returned is eaither a two step or three step progrem
    # This depends on the tm and length of the primers in the
    # original instructions from finnzyme. These do not seem to be

    # Ta calculation for enzymes with dsDNA binding domains like phusion or Pfu-Sso7d
    # https://www.finnzymes.fi/tm_determination.html

    tmf = tm(amplicon.forward_primer.footprint)
    tmr = tm(amplicon.reverse_primer.footprint)

    if tmf >= 69.0 and tmr >= 69.0:
        f = textwrap.dedent(
            r"""
                              |98°C|98°C      |    |tmf:{tmf:.1f}
                              |____|____      |    |tmr:{tmr:.1f}
                              |30s |10s \ 72°C|72°C|{rate}s/kb
                              |    |     \____|____|GC {GC_prod}%
                              |    |     {0:2}:{1:0>2}|5min|{size}bp
                              """.format(
                rate=PfuSso7d_extension_rate,
                tmf=tmf,
                tmr=tmr,
                GC_prod=int(amplicon.gc() * 100),
                size=len(amplicon.seq),
                *map(int, divmod(extension_time_PfuSso7d, 60)),
            )
        ).strip()
    else:
        f = textwrap.dedent(
            r"""
             |98°C|98°C               |    |tmf:{tmf:.1f}
             |____|_____          72°C|72°C|tmr:{tmr:.1f}
             |30s |10s  \ {ta:.1f}°C _____|____|{rate}s/kb
             |    |      \______/{0:2}:{1:0>2}|5min|GC {GC}%
             |    |       10s         |    |{size}bp
             """.format(
                rate=PfuSso7d_extension_rate,
                size=len(amplicon.seq),
                ta=round(
                    ta(
                        amplicon.forward_primer.footprint,
                        amplicon.reverse_primer.footprint,
                        amplicon.seq,
                    ),
                    1,
                ),
                tmf=tmf,
                tmr=tmr,
                GC=int(amplicon.gc() * 100),
                *map(int, divmod(extension_time_PfuSso7d, 60)),
            )
        ).strip()

    return ps(f)


pfu_sso7d_program = dbd_program


def Q5(primer: str, *args, **kwargs):
    """For Q5 Ta they take the lower of the two Tms and add 1C
    (up to 72C). For Phusion they take the lower of the two
    and add 3C (up to 72C).
    """
    raise NotImplementedError


def tmbresluc(primer: str, *args, primerc=500.0, saltc=50, **kwargs):
    """Returns the tm for a primer using a formula adapted to polymerases
    with a DNA binding domain, such as the Phusion polymerase.

    Parameters
    ----------

    primer : string
        primer sequence 5'-3'

    primerc : float
       primer concentration in nM), set to 500.0 nm by default.

    saltc : float, optional
       Monovalent cation concentration in mM, set to 50.0 mM by default.

    thermodynamics : bool, optional
        prints details of the thermodynamic data to stdout. For
        debugging only.

    Returns
    -------
    tm : float
        the tm of the primer

    """

    from . import _thermodynamic_data

    saltc = float(saltc) / 1000
    pri = primerc / 1e9
    dS = -12.4
    dH = -3400

    STR = primer.lower()

    for i in range(len(STR) - 1):
        n1 = ord(STR[i])
        n2 = ord(STR[i + 1])
        dH += _thermodynamic_data.dHBr[n1 - 97][n2 - 97]
        dS += _thermodynamic_data.dSBr[n1 - 97][n2 - 97]

    tm = (
        dH / (1.9872 * math.log(pri / 1600) + dS)
        + (16.6 * math.log(saltc)) / math.log(10)
    ) - 273.15

    return tm


def tm_neb(primer, conc=0.5, prodcode="q5-0"):
    import requests

    """Calculates a single primers melting temp from NEB.

    Parameters
    ----------
    primer1 : str
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    tm : int
        primer melting temperature

    """

    url = "https://tmapi.neb.com/tm"

    params = {"seq1": primer, "conc": conc, "prodcode": prodcode}

    # Mimic a browser request with this headers
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding": "gzip, deflate, br",
        "Accept-Language": "en-US,en;q=0.5",
        "Connection": "keep-alive",
        "Referer": "https://example.com",
    }

    try:
        res = requests.get(url, params=params, headers=headers)
    except requests.exceptions.ConnectionError as e:
        raise requests.exceptions.ConnectionError(
            "Could not connect to NEB API."
        ) from e
    if res.status_code != 200:
        if "error" in res.json():
            raise requests.exceptions.HTTPError(res.status_code, res.json()["error"])
        else:
            raise requests.exceptions.HTTPError(
                res.status_code, res.text
            )  # pragma: no cover
    r = res.json()
    if r["success"]:
        return r["data"]["tm1"]
    else:
        raise requests.exceptions.HTTPError(r["error"])
