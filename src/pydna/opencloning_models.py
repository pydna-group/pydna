# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Optional, Union, Any, ClassVar, Type
from pydantic_core import core_schema

from pydantic import BaseModel, ConfigDict, Field, field_validator

from opencloning_linkml.datamodel import (
    CloningStrategy as _BaseCloningStrategy,
    Primer as PrimerModel,
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
)
from Bio.SeqFeature import Location, LocationParserError
from Bio.Restriction.Restriction import AbstractCut
from typing import List

from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location

from pydna.types import CutSiteType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pydna.dseqrecord import Dseqrecord


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
            id=id(dseqr),
            sequence_file_format="genbank",
            overhang_crick_3prime=dseqr.seq.ovhg,
            overhang_watson_3prime=dseqr.seq.watson_ovhg(),
            file_content=dseqr.format("genbank"),
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
        return _SourceInput(sequence=id(self.sequence))


class AssemblyFragment(SourceInput):
    """
    Represents a fragment in an assembly
    """

    left_location: Optional[Location] = Field(default=None)
    right_location: Optional[Location] = Field(default=None)
    reverse_complemented: bool = Field(
        default=...,
    )

    def to_pydantic_model(self) -> _AssemblyFragment:
        return _AssemblyFragment(
            sequence=id(self.sequence),
            left_location=SequenceLocationStr.from_biopython_location(
                self.left_location
            ),
            right_location=SequenceLocationStr.from_biopython_location(
                self.right_location
            ),
            reverse_complemented=self.reverse_complemented,
        )


class Source(ConfiguredBaseModel):
    input: Optional[list[Union[SourceInput, AssemblyFragment]]] = Field(default=None)

    def to_dict_without_type(self) -> dict:
        return {k: v for k, v in self.model_dump().items() if k != "type"}

    def input_models(self):
        if not self.input:
            return []
        return [fragment.to_pydantic_model() for fragment in self.input]


class AssemblySource(Source):
    circular: bool

    TARGET_MODEL: ClassVar[Type[_AssemblySource]] = _AssemblySource

    def _base_kwargs(self, seq_id: int) -> dict:
        return {
            "id": seq_id,
            "input": self.input_models(),
            "circular": self.circular,
        }

    def extra_kwargs(self) -> dict:
        return {}

    def to_pydantic_model(self, seq_id: int) -> _AssemblySource:
        kwargs = self._base_kwargs(seq_id)
        kwargs.update(self.extra_kwargs())
        return self.TARGET_MODEL(**kwargs)


class RestrictionAndLigationSource(AssemblySource):
    restriction_enzymes: list[AbstractCut]

    TARGET_MODEL: ClassVar[Type[_RestrictionAndLigationSource]] = (
        _RestrictionAndLigationSource
    )

    def extra_kwargs(self) -> dict:
        return {
            "restriction_enzymes": [str(enzyme) for enzyme in self.restriction_enzymes]
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


def cutsite_to_pydantic_model(
    cut_site: CutSiteType | None,
) -> _SequenceCut | _RestrictionSequenceCut | None:
    if cut_site is None:
        return None
    elif isinstance(cut_site[1], AbstractCut):
        return _RestrictionSequenceCut(
            cut_watson=cut_site[0][0],
            overhang=cut_site[0][1],
            restriction_enzyme=str(cut_site[1]),
        )
    else:
        return _SequenceCut(
            cut_watson=cut_site[0][0],
            overhang=cut_site[0][1],
        )


class SequenceCutSource(Source):
    left_edge: CutSiteType | None
    right_edge: CutSiteType | None

    @classmethod
    def from_parent(
        cls, parent: "Dseqrecord", left_edge: CutSiteType, right_edge: CutSiteType
    ):
        return cls(
            input=[SourceInput(sequence=parent)],
            left_edge=left_edge,
            right_edge=right_edge,
        )

    def to_pydantic_model(self, seq_id) -> _SequenceCutSource:
        left_has_enzyme = self.left_edge is not None and isinstance(
            self.left_edge[1], AbstractCut
        )
        right_has_enzyme = self.right_edge is not None and isinstance(
            self.right_edge[1], AbstractCut
        )
        pydantic_class = (
            _RestrictionEnzymeDigestionSource
            if left_has_enzyme or right_has_enzyme
            else _SequenceCutSource
        )
        return pydantic_class(
            id=seq_id,
            input=self.input_models(),
            left_edge=cutsite_to_pydantic_model(self.left_edge),
            right_edge=cutsite_to_pydantic_model(self.right_edge),
        )


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

    def add_primer(self, primer: PrimerModel):
        if primer in self.primers:
            return
        primer.id = self.next_id()
        self.primers.append(primer)

    def next_id(self):
        return (
            max([s.id for s in self.sources + self.sequences + self.primers], default=0)
            + 1
        )

    def add_source_and_sequence(self, source: _Source, sequence: TextFileSequence):
        if source in self.sources:
            if sequence not in self.sequences:
                raise ValueError(
                    (
                        f"Source {source.id} already exists in the cloning strategy, "
                        f"but sequence {sequence.id} it's not its output."
                    )
                )
            return
        new_id = self.next_id()
        source.id = new_id
        self.sources.append(source)
        sequence.id = new_id
        self.sequences.append(sequence)

    def all_children_source_ids(
        self, source_id: int, source_children: list | None = None
    ) -> list[int]:
        """Returns the ids of all source children ids of a source"""
        source = next(s for s in self.sources if s.id == source_id)
        if source_children is None:
            source_children = []

        sources_that_take_output_as_input = [
            s for s in self.sources if source.id in [inp.sequence for inp in s.input]
        ]
        new_source_ids = [s.id for s in sources_that_take_output_as_input]

        source_children.extend(new_source_ids)
        for new_source_id in new_source_ids:
            self.all_children_source_ids(new_source_id, source_children)
        return source_children

    def add_dseqrecord(self, dseqr: "Dseqrecord"):
        existing_ids = {seq.id for seq in self.sequences}
        if id(dseqr) in existing_ids:
            return
        self.sequences.append(TextFileSequence.from_dseqrecord(dseqr))
        if dseqr.source is not None:
            self.sources.append(dseqr.source.to_pydantic_model(id(dseqr)))
            this_source: Source = dseqr.source
            for source_input in this_source.input:
                self.add_dseqrecord(source_input.sequence)
        else:
            self.sources.append(
                _ManuallyTypedSource(id=id(dseqr), input=[], user_input="A")
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
        for source in self.sources:
            source.id = id_mappings[source.id]
            for assembly_fragment in source.input:
                assembly_fragment.sequence = id_mappings[assembly_fragment.sequence]
        for primer in self.primers:
            primer.id = id_mappings[primer.id]

    @classmethod
    def from_dseqrecord(cls, dseqr: "Dseqrecord", description: str = ""):
        cloning_strategy = cls(sources=[], sequences=[], description=description)
        cloning_strategy.add_dseqrecord(dseqr)
        return cloning_strategy

    def model_dump_json(self, *args, **kwargs):
        self.reassign_ids()
        return super().model_dump_json(*args, **kwargs)
