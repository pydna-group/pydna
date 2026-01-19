#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

# PMID: 22645496

# https://www.genome.jp/kegg/catalog/org_list.html

from typing import Dict

_sce_weights = {
    "TTT": 1.0,
    "TTG": 1.0,
    "TAT": 1.0,
    "TAG": 0.479,
    "TGT": 1.0,
    "TGG": 1.0,
    "TCT": 1.0,
    "TCG": 0.364,
    "ATT": 1.0,
    "ATG": 1.0,
    "AAT": 1.0,
    "AAG": 0.736,
    "AGT": 0.602,
    "AGG": 0.433,
    "ACT": 1.0,
    "ACG": 0.393,
    "GTT": 1.0,
    "GTG": 0.488,
    "GAT": 1.0,
    "GAG": 0.422,
    "GGT": 1.0,
    "GGG": 0.252,
    "GCT": 1.0,
    "GCG": 0.292,
    "CTT": 0.451,
    "CTG": 0.386,
    "CAT": 1.0,
    "CAG": 0.444,
    "CGT": 0.3,
    "CGG": 0.082,
    "CCT": 0.738,
    "CCG": 0.289,
    "TTA": 0.962,
    "TTC": 0.706,
    "TAA": 1.0,
    "TAC": 0.787,
    "TGA": 0.643,
    "TGC": 0.588,
    "TCA": 0.795,
    "TCC": 0.605,
    "ATA": 0.59,
    "ATC": 0.57,
    "AAA": 1.0,
    "AAC": 0.696,
    "AGA": 1.0,
    "AGC": 0.415,
    "ACA": 0.876,
    "ACC": 0.628,
    "GTA": 0.533,
    "GTC": 0.533,
    "GAA": 1.0,
    "GAC": 0.538,
    "GGA": 0.456,
    "GGC": 0.409,
    "GCA": 0.765,
    "GCC": 0.595,
    "CTA": 0.493,
    "CTC": 0.2,
    "CAA": 1.0,
    "CAC": 0.571,
    "CGA": 0.141,
    "CGC": 0.122,
    "CCA": 1.0,
    "CCC": 0.37,
}

weights = {"sce": _sce_weights}

# PMID: 6390186
# PMID: 11589713


start: Dict[str, Dict[str, float]] = {
    "sce": {"ATG": 1.000, "TTG": 0.069, "ATA": 0.005},
    "eco": {},
}


# Zhang, S. P., Zubay, G., & Goldman, E. (1991).
# Low-usage codons in Escherichia
# coli, yeast, fruit fly and primates. Gene, 105(1), 61–72.
# https://www.embl.de/pepcore/pepcore_services/cloning/choice_expression_systems/codons8
# AGG ACG         CGA CGG CGC CCG CTC GCG


rare_codons = {
    "sce": ["CGA", "CGG", "CGC", "CCG", "CTC", "GCG"],
    "eco": ["AGG", "AGA", "ATA", "CTA", "CGA", "CGG", "CCC", "TCG"],
}

stop: Dict[str, Dict[str, float]] = {
    "sce": {"TAA": 0.470, "TAG": 0.230, "TGA": 0.300},
    "eco": {},
}


n_end: Dict[str, Dict[str, str]] = {
    "sce": {
        "Val": ">30 h",
        "Met": ">30 h",
        "Gly": ">30 h",
        "Pro": ">5 h",
        "Ala": ">30 h",
        "Ser": ">30 h",
        "Thr": ">30 h",
        "Cys": ">30 h",
        "Ile": "30 min",
        "Glu": "30 min",
        "His": "3 min",
        "Tyr": "10 min",
        "Gln": "10 min",
        "Asp": "3 min",
        "Asn": "3 min",
        "Phe": "3 min",
        "Leu": "3 min",
        "Trp": "3 min",
        "Lys": "3 min",
        "Arg": "2 min",
    },
    "eco": {},
}
