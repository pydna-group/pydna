# -*- coding: utf-8 -*-
"""Homologous recombination between shared sequence.

A genome is either opened to take up an insert (*integration*) or resolved at
two of its own homologous stretches (*excision or inversion*). In-vivo assembly,
which relies on the same chemistry between separate fragments, is declared in
:mod:`pydna.methods.ligation`.
"""

from __future__ import annotations

import warnings

from pydna.assembly import common_handle_insertion_fragments, common_sub_strings
from pydna.core.dseqrecord import Dseqrecord
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import HomologousRecombinationSource

__all__ = [
    "homologous_recombination_integration",
    "homologous_recombination_excision_or_inversion",
    "homologous_recombination_excision",
]


def _genome_and_inserts(params: dict) -> list[Dseqrecord]:
    """The genome first, then the inserts, as the integration shape expects."""
    return common_handle_insertion_fragments(
        params.pop("genome"), params.pop("inserts")
    )


def _genome_only(params: dict) -> list[Dseqrecord]:
    return [params.pop("genome")]


# --- homologous recombination ------------------------------------------------

homologous_recombination_integration = register(
    Method(
        name="homologous_recombination_integration",
        doc="""Returns the products resulting from the integration of an insert (or inserts joined
    through in vivo recombination) into the genome through homologous recombination.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord]
        DNA fragment(s) to insert
    limit : int, optional
        Minimum homology length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules


    Examples
    --------

    Below an example with a single insert.

    >>> from pydna.methods import homologous_recombination_integration
    >>> from pydna.core.dseqrecord import Dseqrecord
    >>> homology = "AAGTCCGTTCGTTTTACCTG"
    >>> genome = Dseqrecord(f"aaaaaa{homology}ccccc{homology}aaaaaa")
    >>> insert = Dseqrecord(f"{homology}gggg{homology}")
    >>> products = homologous_recombination_integration(genome, [insert], 20)
    >>> str(products[0].seq)
    'aaaaaaAAGTCCGTTCGTTTTACCTGggggAAGTCCGTTCGTTTTACCTGaaaaaa'

    Below an example with two inserts joined through homology.

    >>> homology2 = "ATTACAGCATGGGAAGAAAGA"
    >>> insert_1 = Dseqrecord(f"{homology}gggg{homology2}")
    >>> insert_2 = Dseqrecord(f"{homology2}cccc{homology}")
    >>> products = homologous_recombination_integration(genome, [insert_1, insert_2], 20)
    >>> str(products[0].seq)
    'aaaaaaAAGTCCGTTCGTTTTACCTGggggATTACAGCATGGGAAGAAAGAccccAAGTCCGTTCGTTTTACCTGaaaaaa'
        """,
        shape=Shape.INTEGRATION,
        algorithm=common_sub_strings,
        source=HomologousRecombinationSource,
        limit=40,
        prepare_inputs=_genome_and_inserts,
        positional=("genome", "inserts", "limit"),
        summary="Integrate insert(s) into a genome by homologous recombination.",
    )
)

homologous_recombination_excision_or_inversion = register(
    Method(
        name="homologous_recombination_excision_or_inversion",
        doc="""Returns the products resulting from the excision of a fragment from the genome through
    homologous recombination.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    limit : int, optional
        Minimum homology length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List containing excised plasmid and remaining genome sequence

    Examples
    --------

    Example of a homologous recombination event, where a plasmid is excised from the
    genome (circular sequence of 25 bp), and that part is removed from the genome,
    leaving a shorter linear sequence (32 bp).

    >>> from pydna.methods import homologous_recombination_excision
    >>> from pydna.core.dseqrecord import Dseqrecord
    >>> homology = "AAGTCCGTTCGTTTTACCTG"
    >>> genome = Dseqrecord(f"aaaaaa{homology}ccccc{homology}aaaaaa")
    >>> products = homologous_recombination_excision(genome, 20)
    >>> products
    [Dseqrecord(o25), Dseqrecord(-32)]
        """,
        shape=Shape.EXCISION,
        algorithm=common_sub_strings,
        source=HomologousRecombinationSource,
        limit=40,
        prepare_inputs=_genome_only,
        positional=("genome", "limit"),
        summary="Excise or invert a genome region between homologous sequences.",
    )
)


def homologous_recombination_excision(
    genome: Dseqrecord, limit: int = 40
) -> list[Dseqrecord]:
    """Deprecated alias of homologous_recombination_excision_or_inversion."""
    warnings.warn(
        "`homologous_recombination_excision` is deprecated and will be removed in a future "
        "version; use `homologous_recombination_excision_or_inversion` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return homologous_recombination_excision_or_inversion(genome, limit)
