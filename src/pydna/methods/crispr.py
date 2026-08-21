#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

"""
Utilities for CRISPR/Cas target searching and protospacer extraction.

"""

import re
from typing import Type
from typing import TYPE_CHECKING
from typing import List
from typing import TypeVar
import warnings
from pydna.assembly import common_handle_insertion_fragments, common_sub_strings
from pydna.core.dseqrecord import Dseqrecord
from pydna.core.nuclease import _cas
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import CRISPRSource, SourceInput
from pydna.utils import create_location, location_boundaries

if TYPE_CHECKING:  # pragma: no cover
    from pydna.core.dseqrecord import Dseqrecord

DseqrecordType = TypeVar("DseqrecordType", bound="Dseqrecord")


class cas9(_cas):
    """docstring.

    .. code-block::


            fst5              --|fst3
            |----------------
                                 PAM
       5'-NNGGAAGAGTAATACACTA-AAANGGNN-3'
          ||||||||||||||||||| ||||||||
       3'-NNCCTTCTCATTATGTGAT-TTTNCCNN-5'
            ||||||||||||||||| |||
         5'-GGAAGAGTAATACACTA AAA-g-u-a-a-g-g-3'  Scaffold (lower case)
            ---gRNA spacer------- u-a
                                  u-a
                                  u-a
                                  u-a
                                  a-u
                                  g-u-g
                                  a    a
                                  g-c-a
                                  c-g
                                  u-a
                                  a-u
                                  g   a
                                   a-a
    """

    scaffold: str = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    pam: str = ".GG"
    size: int = 20
    fst5: int = 17
    fst3: int = -3
    ovhg: int = fst5 - (size + fst3)

    def search(self, dna, linear: bool = True) -> List[int]:
        """
        Search for Cas9 target sites in a DNA sequence.

        Args:
            dna: string, Bio.Seq.Seq or pydna.core.dseq.Dseq
            linear: Whether the DNA is linear or circular.

        Returns:
            A list of cut site positions.
        """
        from pydna.core.dseqrecord import Dseqrecord
        from pydna.sequence_regex import dseqrecord_finditer

        if not hasattr(dna, "_data"):
            raise TypeError
        results: List[int] = []
        query = Dseqrecord(dna, circular=(not linear))
        matches_fwd = dseqrecord_finditer(self.compsite, query)
        matches_rev = dseqrecord_finditer(self.compsite, query.reverse_complement())
        for mobj in matches_fwd:
            results.append((mobj.start() + self.fst5 + 1) % len(dna))
        for mobj in matches_rev:
            results.append((len(dna) - (mobj.start() + self.fst5) + 1) % len(dna))
        return results


def protospacer(guide_construct: DseqrecordType, cas: Type[_cas] = cas9) -> List[str]:
    """
    Extract protospacer sequences from a guide construct. This can for example
    be a plasmid containing the guide construct. This function returns a
    list since a several protospacers can be present.

    Args:
        guide_construct: Sequence construct containing protospacer and scaffold.
        cas: CRISPR nuclease class defining spacer size and scaffold.

    Returns:
        A list of protospacer sequences found in Watson and Crick orientations.
    """
    if guide_construct.circular:
        total_length = cas.size + len(cas.scaffold)
        guide_construct = guide_construct[:] + guide_construct[: total_length - 1]

    result = []

    for s in guide_construct.seq.watson.upper(), guide_construct.seq.crick.upper():
        result.extend(
            mobj.group("ps")
            for mobj in re.finditer(
                f"(?P<ps>.{{{cas.size}}})(?:{cas.scaffold})",
                s,
            )
        )

    return result


def _validate(inputs, params):
    if len(params.get("guides") or []) == 0:
        raise ValueError("At least one guide RNA is required for CRISPR integration")


def _prepare_inputs(params: dict) -> list[Dseqrecord]:
    return common_handle_insertion_fragments(
        params.pop("genome"), params.pop("inserts")
    )


def _guides_cut_inside_repair(products, inputs, params):
    """Keep products whose repaired region contains every guide's cut site.

    Each surviving product records the guides that made it, which is what lets
    the strategy be replayed.
    """
    genome = inputs[0]
    guides = params["guides"]

    # First we collect the positions where the guides cut
    guide_cuts = []
    for guide in guides:
        enzyme = cas9(str(guide.seq))
        possible_cuts = genome.seq.get_cutsites(enzyme)
        if len(possible_cuts) == 0:
            raise ValueError(
                f"Could not find Cas9 cutsite in the target sequence using the guide: {guide.name}"
            )
        # Keep only the position of the cut
        possible_cuts = [cut[0] for (cut, _) in possible_cuts]
        guide_cuts.append(possible_cuts)

    valid_products = []
    for product in products:
        # The second element of product.source.input is conventionally the insert/repair
        # fragment; the other two (first and third) are the two bits of the genome.
        repair_start = location_boundaries(product.source.input[0].right_location)[0]
        # +1 because the position of the cut marks the boundary
        repair_end = location_boundaries(product.source.input[2].left_location)[1] + 1
        repair_location = create_location(repair_start, repair_end, len(genome))

        some_cuts_inside_repair = []
        all_cuts_inside_repair = []
        for cut_group in guide_cuts:
            cuts_in_repair = [cut for cut in cut_group if cut in repair_location]
            some_cuts_inside_repair.append(len(cuts_in_repair) != 0)
            all_cuts_inside_repair.append(len(cuts_in_repair) == len(cut_group))

        if all(some_cuts_inside_repair):
            used_guides = [g for i, g in enumerate(guides) if all_cuts_inside_repair[i]]
            # Record the guides that produced this product
            product.source.input.extend([SourceInput(sequence=g) for g in used_guides])
            valid_products.append(product)
            if not all(all_cuts_inside_repair):
                raise ValueError(
                    "Some guides cut outside the repair region, please check the guides"
                )

    if len(valid_products) != len(products):
        warnings.warn(
            "Some recombination products were discarded because they had off-target cuts",
            category=UserWarning,
            stacklevel=2,
        )
    return valid_products


crispr_integration = register(
    Method(
        name="crispr_integration",
        doc="""
    Returns the products for CRISPR integration.

    Parameters
    ----------
    genome : Dseqrecord
        Target genome sequence
    inserts : list[Dseqrecord]
        DNA fragment(s) to insert
    guides : list[Primer]
        List of guide RNAs as Primer objects. This may change in the future.
    limit : int, optional
        Minimum overlap length required, by default 40

    Returns
    -------
    list[Dseqrecord]
        List of integrated DNA molecules

    Examples
    --------

    >>> from pydna.core.dseqrecord import Dseqrecord
    >>> from pydna.methods import crispr_integration
    >>> from pydna.primer import Primer
    >>> genome = Dseqrecord("aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac", name="genome")
    >>> insert = Dseqrecord("aaccggttAAAAAAAAAttcaaagcac", name="insert")
    >>> guide = Primer("ttcaatgcaaacagtaatga", name="guide")
    >>> product, *_ = crispr_integration(genome, [insert], [guide], 8)
    >>> product
    Dseqrecord(-27)
        """,
        shape=Shape.INTEGRATION,
        algorithm=common_sub_strings,
        source=CRISPRSource,
        limit=40,
        validate=_validate,
        filter_products=_guides_cut_inside_repair,
        prepare_inputs=_prepare_inputs,
        positional=("genome", "inserts", "guides", "limit"),
        summary="Integrate insert(s) at a Cas9 cut directed by guide RNAs.",
    )
)
