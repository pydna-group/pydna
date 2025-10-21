# -*- coding: utf-8 -*-
"""
This module provides classes that roughly map to the `OpenCloning <https://opencloning.org>`_
data model, which is defined using `LinkML <https://linkml.io>`, and available as a python
package `opencloning-linkml <https://pypi.org/project/opencloning-linkml/>`_. These classes
are documented there, and the ones in this module essentially replace the fields pointing to
sequences and primers (which use ids in the data model) to ``Dseqrecord`` and ``Primer``
objects, respectively. Similarly, it uses Location from ``Biopython`` instead of a string,
which is what the data model uses.

When using pydna to plan cloning, it stores the provenance of ``Dseqrecord`` objects in
their ``source`` attribute. Not all methods generate sources so far, so refer to the
documentation notebooks for examples on how to use this feature. The ``history`` method of
``Dseqrecord`` objects can be used to get a string representation of the provenance of the
sequence. You can also use the ``CloningStrategy`` class to create a JSON representation of
the cloning strategy. That ``CloningStrategy`` can be loaded in the OpenCloning web interface
to see a representation of the cloning strategy.

"""
from __future__ import annotations

from typing import Optional, Union, Any, ClassVar, Type
from pydantic_core import core_schema
from contextlib import contextmanager
from threading import local
import re

from pydantic import BaseModel, ConfigDict, Field, field_validator

from opencloning_linkml.datamodel import (
    CloningStrategy as _BaseCloningStrategy,
    Primer as _PrimerModel,
    Source as _Source,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    SourceInput as _SourceInput,
    AssemblyFragment as _AssemblyFragment,
    ManuallyTypedSource as _ManuallyTypedSource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    SequenceCutSource as _SequenceCutSource,
    RestrictionSequenceCut as _RestrictionSequenceCut,
    SequenceCut as _SequenceCut,
    InFusionSource as _InFusionSource,
    OverlapExtensionPCRLigationSource as _OverlapExtensionPCRLigationSource,
    InVivoAssemblySource as _InVivoAssemblySource,
    LigationSource as _LigationSource,
    GatewaySource as _GatewaySource,
    GatewayReactionType,
    HomologousRecombinationSource as _HomologousRecombinationSource,
    CreLoxRecombinationSource as _CreLoxRecombinationSource,
    PCRSource as _PCRSource,
    CRISPRSource as _CRISPRSource,
)
from Bio.SeqFeature import Location, LocationParserError
from Bio.Restriction.Restriction import AbstractCut
import networkx as nx
from typing import List

from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location

from pydna.types import (
    CutSiteType,
    SubFragmentRepresentationAssembly,
    EdgeRepresentationAssembly,
)
from pydna.utils import create_location
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer


# Thread-local storage for ID strategy
_thread_local = local()


@contextmanager
def id_mode(use_python_internal_id: bool = True):
    """Context manager that is used to determine how ids are assigned to objects when
    mapping them to the OpenCloning data model. If ``use_python_internal_id`` is True,
    the built-in python ``id()`` function is used to assign ids to objects. That function
    produces a unique integer for each object in python, so it's guaranteed to be unique.
    If ``use_python_internal_id`` is False, the object's ``.id`` attribute (must be a string integer)
    is used to assign ids to objects. This is useful when the objects already have meaningful ids,
    and you want to keep references to them in ``SourceInput`` objects (which sequences and
    primers are used in a particular source).

    Parameters
    ----------
    use_python_internal_id: bool
        If True, use Python's built-in id() function.
        If False, use the object's .id attribute (must be a string integer).

    Examples
    --------
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.opencloning_models import get_id, id_mode
    >>> dseqr = Dseqrecord("ATGC")
    >>> dseqr.name = "my_sequence"
    >>> dseqr.id = "123"
    >>> get_id(dseqr) == id(dseqr)
    True
    >>> with id_mode(use_python_internal_id=False):
    ...     get_id(dseqr)
    123
    """
    old_value = getattr(_thread_local, "use_python_internal_id", True)
    _thread_local.use_python_internal_id = use_python_internal_id
    try:
        yield
    finally:
        _thread_local.use_python_internal_id = old_value


def get_id(obj: "Primer" | "Dseqrecord") -> int:
    """Get ID using the current strategy from thread-local storage (see id_mode)
    Parameters
    ----------
    obj: Primer | Dseqrecord
        The object to get the id of

    Returns
    -------
    int: The id of the object

    """
    use_python_internal_id = getattr(_thread_local, "use_python_internal_id", True)
    if use_python_internal_id:
        return id(obj)
    if not isinstance(obj.id, str) or not obj.id.isdigit():
        raise ValueError(
            f"If use_python_internal_id is False, id must be a string representing an integer, "
            f"but object {obj} has an invalid id: {obj.id}"
        )
    return int(obj.id)


class SequenceLocationStr(str):
    """A string representation of a sequence location, genbank-like."""

    # TODO: this should handle origin-spanning simple locations (splitted)
    @classmethod
    def from_biopython_location(cls, location: Location):
        return cls(format_feature_location(location, None))

    def to_biopython_location(self) -> Location:
        return Location.fromstring(self)

    @classmethod
    def field_validator(cls, v):
        if isinstance(v, str):
            value = cls(v)
            try:
                value.to_biopython_location()
            except LocationParserError as err:
                raise ValueError(f"Location {v!r} is not a valid location") from err
            return value
        raise ValueError(f"Location must be a string or a {cls.__name__}")

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        source_type,
        handler,
    ) -> core_schema.CoreSchema:
        """Generate Pydantic core schema for SequenceLocationStr."""
        return core_schema.with_info_after_validator_function(
            cls._validate,
            core_schema.str_schema(),
        )

    @classmethod
    def _validate(cls, value: str, info):
        """Validate and create SequenceLocationStr instance."""
        return cls.field_validator(value)

    @classmethod
    def from_start_and_end(
        cls, start: int, end: int, seq_len: int | None = None, strand: int | None = 1
    ):
        return cls.from_biopython_location(create_location(start, end, seq_len, strand))


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        validate_assignment=True,
        validate_default=True,
        extra="forbid",
        arbitrary_types_allowed=True,
        use_enum_values=True,
        strict=False,
    )
    pass


class TextFileSequence(_TextFileSequence):

    @classmethod
    def from_dseqrecord(cls, dseqr: "Dseqrecord"):
        return cls(
            id=get_id(dseqr),
            sequence_file_format="genbank",
            overhang_crick_3prime=dseqr.seq.ovhg,
            overhang_watson_3prime=dseqr.seq.watson_ovhg(),
            file_content=dseqr.format("genbank"),
        )


class PrimerModel(_PrimerModel):

    @classmethod
    def from_primer(cls, primer: "Primer"):
        return cls(
            id=get_id(primer),
            name=primer.name,
            sequence=str(primer.seq),
        )


class SourceInput(ConfiguredBaseModel):
    sequence: object

    @field_validator("sequence")
    @classmethod
    def _validate_sequence_field(cls, value: Any):
        """Separate validation to avoid circular imports."""

        from pydna.dseqrecord import Dseqrecord
        from pydna.primer import Primer

        if isinstance(value, (Dseqrecord, Primer)):
            return value
        module = type(value).__module__
        name = type(value).__name__
        raise TypeError(f"sequence must be Dseqrecord or Primer; got {module}.{name}")

    def to_pydantic_model(self) -> _SourceInput:
        return _SourceInput(sequence=get_id(self.sequence))


class AssemblyFragment(SourceInput):

    left_location: Optional[Location] = Field(default=None)
    right_location: Optional[Location] = Field(default=None)
    reverse_complemented: bool

    @staticmethod
    def from_biopython_location(location: Location | None):
        if location is None:
            return None
        return SequenceLocationStr.from_biopython_location(location)

    def to_pydantic_model(self) -> _AssemblyFragment:
        return _AssemblyFragment(
            sequence=get_id(self.sequence),
            left_location=self.from_biopython_location(self.left_location),
            right_location=self.from_biopython_location(self.right_location),
            reverse_complemented=self.reverse_complemented,
        )


class Source(ConfiguredBaseModel):
    input: list[Union[SourceInput, AssemblyFragment]] = Field(default_factory=list)
    TARGET_MODEL: ClassVar[Type[_Source]] = _Source

    def input_models(self):
        return [fragment.to_pydantic_model() for fragment in self.input]

    def _kwargs(self, seq_id: int) -> dict:
        return {
            "id": seq_id,
            "input": self.input_models(),
        }

    def to_pydantic_model(self, seq_id: int):
        kwargs = self._kwargs(seq_id)
        return self.TARGET_MODEL(**kwargs)

    def add_to_history_graph(self, history_graph: nx.DiGraph, seq: "Dseqrecord"):
        """
        Add the source to the history graph.

        It does not use the get_id function, because it just uses it to have unique identifiers
        for graph nodes, not to store them anywhere.
        """
        from pydna.dseqrecord import Dseqrecord

        history_graph.add_node(id(seq), label=f"{seq.name} ({repr(seq)})")
        history_graph.add_node(id(self), label=str(self.TARGET_MODEL.__name__))
        history_graph.add_edge(id(seq), id(self))
        for fragment in self.input:
            fragment_seq = fragment.sequence
            # This could be a Primer as well, which doesn't have a source
            if isinstance(fragment_seq, Dseqrecord) and fragment_seq.source is not None:
                fragment_seq.source.add_to_history_graph(history_graph, fragment_seq)
            else:
                history_graph.add_node(
                    id(fragment_seq),
                    label=f"{fragment_seq.name} ({repr(fragment_seq)})",
                )
            history_graph.add_edge(id(self), id(fragment_seq))

    def history_string(self, seq: "Dseqrecord"):
        """
        Returns a string representation of the cloning history of the sequence.
        See dseqrecord.history() for examples.
        """
        history_graph = nx.DiGraph()
        self.add_to_history_graph(history_graph, seq)
        out_str = "\n".join(
            nx.generate_network_text(history_graph, with_labels=True, sources=[id(seq)])
        )
        # Replace special characters with normal ASCII characters, otherwise
        # notebooks are not properly rendered in GitHub
        replace_dict = {
            "╙": "*",  # root node
            "└": "+",  # corner
            "├": "|",  # junction
            "╼": ">",  # arrow
            "─": "-",  # horizontal line
            "╾": "~",  # junction
            "│": "|",  # vertical line
        }
        out_str = re.sub(
            rf"{'|'.join(replace_dict.keys())}",
            lambda m: replace_dict[m.group()],
            out_str,
        )
        return out_str


class AssemblySource(Source):
    circular: bool

    TARGET_MODEL: ClassVar[Type[_AssemblySource]] = _AssemblySource

    def _kwargs(self, seq_id: int) -> dict:
        return {
            **super()._kwargs(seq_id),
            "circular": self.circular,
        }

    def to_pydantic_model(self, seq_id: int):
        return self.TARGET_MODEL(**self._kwargs(seq_id))

    @classmethod
    def from_subfragment_representation(
        cls,
        assembly: SubFragmentRepresentationAssembly,
        fragments: list["Dseqrecord"],
        is_circular: bool,
    ):

        input_list = []
        for f_index, loc1, loc2 in assembly:
            input_list.append(
                AssemblyFragment(
                    sequence=fragments[abs(f_index) - 1],
                    left_location=loc1,
                    right_location=loc2,
                    reverse_complemented=f_index < 0,
                )
            )

        return AssemblySource(input=input_list, circular=is_circular)

    def from_edge_representation(
        cls,
        edge_representation: EdgeRepresentationAssembly,
        fragments: list["Dseqrecord"],
        is_circular: bool,
    ):
        pass


class RestrictionAndLigationSource(AssemblySource):
    restriction_enzymes: list[AbstractCut]

    TARGET_MODEL: ClassVar[Type[_RestrictionAndLigationSource]] = (
        _RestrictionAndLigationSource
    )

    def _kwargs(self, seq_id: int) -> dict:
        return {
            **super()._kwargs(seq_id),
            "restriction_enzymes": [str(enzyme) for enzyme in self.restriction_enzymes],
        }


class GibsonAssemblySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_GibsonAssemblySource]] = _GibsonAssemblySource


class InFusionSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_InFusionSource]] = _InFusionSource


class OverlapExtensionPCRLigationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_OverlapExtensionPCRLigationSource]] = (
        _OverlapExtensionPCRLigationSource
    )


class InVivoAssemblySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_InVivoAssemblySource]] = _InVivoAssemblySource


class LigationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_LigationSource]] = _LigationSource


class GatewaySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_GatewaySource]] = _GatewaySource
    reaction_type: GatewayReactionType
    greedy: bool = Field(default=False)

    def _kwargs(self, seq_id: int) -> dict:
        return {
            **super()._kwargs(seq_id),
            "reaction_type": self.reaction_type,
            "greedy": self.greedy,
        }


class HomologousRecombinationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_HomologousRecombinationSource]] = (
        _HomologousRecombinationSource
    )


class CRISPRSource(HomologousRecombinationSource):
    TARGET_MODEL: ClassVar[Type[_CRISPRSource]] = _CRISPRSource


class CreLoxRecombinationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_CreLoxRecombinationSource]] = (
        _CreLoxRecombinationSource
    )


class PCRSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_PCRSource]] = _PCRSource
    add_primer_features: bool = Field(default=False)

    def _kwargs(self, seq_id: int) -> dict:
        return {
            **super()._kwargs(seq_id),
            "add_primer_features": self.add_primer_features,
        }


class SequenceCutSource(Source):
    left_edge: CutSiteType | None
    right_edge: CutSiteType | None

    BASE_MODEL: ClassVar[Type[_SequenceCutSource]] = _SequenceCutSource
    ENZYME_MODEL: ClassVar[Type[_RestrictionEnzymeDigestionSource]] = (
        _RestrictionEnzymeDigestionSource
    )

    @staticmethod
    def _cutsite_to_model(cut_site: CutSiteType | None):
        if cut_site is None:
            return None
        watson, overhang = cut_site[0]
        enzyme_or_none = cut_site[1]
        if isinstance(enzyme_or_none, AbstractCut):
            return _RestrictionSequenceCut(
                cut_watson=watson,
                overhang=overhang,
                restriction_enzyme=str(enzyme_or_none),
            )
        return _SequenceCut(cut_watson=watson, overhang=overhang)

    @classmethod
    def from_parent(
        cls, parent: "Dseqrecord", left_edge: CutSiteType, right_edge: CutSiteType
    ):
        return cls(
            input=[SourceInput(sequence=parent)],
            left_edge=left_edge,
            right_edge=right_edge,
        )

    def _has_enzyme(self) -> bool:
        def has_enzyme(edge):
            return edge is not None and isinstance(edge[1], AbstractCut)

        return has_enzyme(self.left_edge) or has_enzyme(self.right_edge)

    def _target_model(self):
        return self.ENZYME_MODEL if self._has_enzyme() else self.BASE_MODEL

    def _kwargs(self, seq_id: int) -> dict:
        return {
            **super()._kwargs(seq_id),
            "left_edge": self._cutsite_to_model(self.left_edge),
            "right_edge": self._cutsite_to_model(self.right_edge),
        }

    def to_pydantic_model(self, seq_id: int):
        return self._target_model()(**self._kwargs(seq_id))


class CloningStrategy(_BaseCloningStrategy):

    # For now, we don't add anything, but the classes will not have the new
    # methods if this is used
    # It will be used for validation for now
    primers: Optional[List[PrimerModel]] = Field(
        default_factory=list,
        description="""The primers that are used in the cloning strategy""",
        json_schema_extra={
            "linkml_meta": {"alias": "primers", "domain_of": ["CloningStrategy"]}
        },
    )

    def add_primer(self, primer: "Primer"):
        existing_ids = {seq.id for seq in self.primers}
        if get_id(primer) in existing_ids:
            return
        self.primers.append(PrimerModel.from_primer(primer))

    def add_dseqrecord(self, dseqr: "Dseqrecord"):
        from pydna.dseqrecord import Dseqrecord

        existing_ids = {seq.id for seq in self.sequences}
        if get_id(dseqr) in existing_ids:
            return
        self.sequences.append(TextFileSequence.from_dseqrecord(dseqr))
        if dseqr.source is not None:
            self.sources.append(dseqr.source.to_pydantic_model(get_id(dseqr)))
            this_source: Source = dseqr.source
            for source_input in this_source.input:
                if isinstance(source_input.sequence, Dseqrecord):
                    self.add_dseqrecord(source_input.sequence)
                else:
                    self.add_primer(source_input.sequence)
        else:
            self.sources.append(
                _ManuallyTypedSource(id=get_id(dseqr), input=[], user_input="A")
            )

    def reassign_ids(self):
        all_ids = (
            {seq.id for seq in self.sequences}
            | {source.id for source in self.sources}
            | {primer.id for primer in self.primers}
        )
        id_mappings = {id: i + 1 for i, id in enumerate(sorted(all_ids))}
        for seq in self.sequences:
            seq.id = id_mappings[seq.id]
        for primer in self.primers:
            primer.id = id_mappings[primer.id]
        for source in self.sources:
            source.id = id_mappings[source.id]
            for assembly_fragment in source.input:
                assembly_fragment.sequence = id_mappings[assembly_fragment.sequence]

    @classmethod
    def from_dseqrecords(cls, dseqrs: list["Dseqrecord"], description: str = ""):
        cloning_strategy = cls(sources=[], sequences=[], description=description)
        for dseqr in dseqrs:
            cloning_strategy.add_dseqrecord(dseqr)
        return cloning_strategy

    def model_dump_json(self, *args, **kwargs):
        if getattr(_thread_local, "use_python_internal_id", True):
            # Make a deep copy of the cloning strategy and reassign ids
            cs = self.__deepcopy__()
            cs.reassign_ids()
            return super(CloningStrategy, cs).model_dump_json(*args, **kwargs)
        return super().model_dump_json(*args, **kwargs)

    def model_dump(self, *args, **kwargs):
        if getattr(_thread_local, "use_python_internal_id", True):
            cs = self.__deepcopy__()
            cs.reassign_ids()
            return super(CloningStrategy, cs).model_dump(*args, **kwargs)
        return super().model_dump(*args, **kwargs)
